/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' mapRead SOF: Start Of File
'   - maps a read to a reference
'   * .h st01: chains_mapRead
'     - holds kmer chains build from kmerChain_mapRead
'      o fun01: blank_chains_mapRead
'        - blanks a chains_mapRead struct
'      o fun02: init_chains_mapRead
'        - initialize a chains_mapRead struct
'      o fun03: freeStack_chains_mapRead
'        - frees all variables in a chains_mapRead struct
'      o fun04: freeHeap_chains_mapRead
'        - frees a chains_mapRead struct
'      o fun05: setup_chains_mapRead
'        - sets up memory for a chains_mapRead struct
'      o fun06: mk_chains_mapRead
'        - makes a heap allocated chains_mapRead struct
'      o fun07: realloc_chains_mapRead
'        - reallocates memroy in a chains_mapRead struct
'   * .h st02: set_mapRead
'     - settings for a aln_mapRead
'     o fun08: blank_set_mapRead
'       - blank a set_mapRead struct
'     o fun09: init_set_mapRead
'       - initializes a set_mapRead struct
'     o fun10: freeStack_set_mapRead
'       - frees variables in a set_mapRead struct
'     o fun11: freeHeap_set_mapRead
'       - frees a set_mapRead struct
'     o fun12: setup_set_mapRead
'       - allocates memory for structs in a set_mapRead
'   * .h st03: ref_mapRead
'     - holds reference variables used in read mapping
'     o fun13: blank_ref_mapRead
'       - blanks (sets to defaults) a ref_mapRead struct
'     o fun14: init_ref_mapRead
'       - initialize a ref_mapRead structure
'     o fun15: freeStack_ref_mapRead
'       - frees variables in a ref_mapRead struct
'     o fun16: freeHeap_ref_mapRead
'       - frees a ref_mapRead struct
'     o fun17: addRef_ref_mapRead
'       - adds a reference sequence to ref_map read struct
'     o fun18: mk_ref_mapRead
'       - makes and initializes a ref_mapRead structure
'     o fun19: getRef_ref_mapRead
'       - gets and stores next reference in ref_mapRead
'     o fun20: getRefBin_ref_mapRead
'       - gets a binary file with reference
'     o fun21: writeRefBin_ref_mapRead
'       - writes a ref_mapRead struct to a binary file
'   * .h st04: aln_mapRead
'     - holds non-reference variables used in read mapping
'     o fun22: blank_aln_mapRead
'       - blanks an aln_mapRead struct
'     o fun23: init_aln_mapRead
'       - initializes an aln_mapRead struct
'     o fun24: freeStack_aln_mapRead
'       - frees variables in an aln_mapRead struct
'     o fun25: freeHeap_aln_mapRead
'       - frees an aln_mapRead struct
'     o fun26: setup_aln_mapRead
'       - setups (memory allocate) an aln_mapRead struct
'     o fun27: mk_aln_mapRead
'       - makes and returns heap allocated aln_mapRead
'   * .c st05: hit_mapRead
'     - private structure to create a linked list of chain
'       index's that I have already aligned
'   o fun28: sortKmerIndex_mapRead
'     - sorts kmer index array by a kmer array
'   o .c fun29: findKmer_mapRead
'     - finds first signed int in a range
'   o .c fun30: findChain_chains_mapRead
'     - finds kmer in range from kmerChain_mapRead return
'   o fun31: kmerChain_mapRead
'     - merge kmers into chains (longer kmers)
'   o .c fun32: numQryChains_mapRead
'     - gets the number of chains assigned to a single
'       query kmer
'   o fun33: mergeChains_mapRead
'     - merges kmer chains from fun31 (kmerChain) into
'       the chain the maximizes base count
'   o fun34: mergeToSam_mapRead
'     - convert output from mergeChains_map read to sam
'       file entry (mapq is set to 0)
'   o fun35: scoreSubAln_mapRead
'     - scores the aligment in samEntry and finds best
'       sub-aligment
'   o fun37: align_mapRead
'     - maps read to reference
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "mapRead.h"

#include <stdio.h>

#include "../genLib/numToStr.h"
#include "../genLib/ulCp.h"
#include "../genLib/shellSort.h"

#include "../genBio/kmerFun.h"
#include "../genBio/seqST.h"
#include "../genBio/samEntry.h"

#include "alnSet.h"
#include "dirMatrix.h"
#include "needle.h"
#include "water.h"

/*.h files only*/
#include "defsMapRead.h" /*default settings*/
#include "../genLib/genMath.h" /*.h max/min macros only*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden libraries
!   - .c  #include "../genLib/base10Str.h"
!   - .c  #include "../genLib/charCp.h"
!   - .c  #include "../genLib/strAry.h"
!   - .c  #include "../genLib/fileFun.h"
!   - .c  #include "../genAln/indexToCoord.h"
!   - .h  #include "../genBio/ntTo5Bit.h"
!   - .h  #include "../genLib/endLine.h"
!   - .h  #include "../genBio/kmerBit.h"
!   - .h  #include "../genAln/alnDefs.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_minSL_mapRead ((signed long) ( (unsigned long) -1 ^ ((unsigned long) -1 >> 1) ) )
   /*Logic: smallest long possible
   `   - maxSL: (unsigned long) -1 >> 1:
   `     o max value for signed long; removes negative bit
   `   - (unsigned long) -1 ^ maxSL:
   `     o clears all bits except the negative (min value)
   */

/*-------------------------------------------------------\
| ST05: hit_mapRead
|   - private structure to create a linked list of chain
|     index's that I have already aligned
\-------------------------------------------------------*/
typedef struct hit_mapRead
{
   signed int indexSI;
   struct hit_mapRead *nextST;
}hit_mapRead;

/*-------------------------------------------------------\
| Fun01: blank_chains_mapRead
|   - blanks a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct to blank
| Output:
|   - Modifies:
|     o refKmerAryUI and qryKmerAryUI in chainsSTPtr to
|       have all values set to -1
|     o chainLenArySI and scoreArySI in chainsSTPtr to
|       have all values set to 0
|     o indexArySI in chainsSTPtr to be filled with -1's
|     o lenChainSI and lenIndexSL in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
blank_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to blank*/
){
   signed long slElm = 0;

   if(! chainsSTPtr)
      return; /*nothing to blank*/

   if(chainsSTPtr->refArySI)
   { /*If: need to blank chain arrays*/

      slElm = 0;

      while(slElm < chainsSTPtr->lenChainSI)
      { /*Loop: blank chain arrays*/
         chainsSTPtr->refArySI[slElm] = -1;
         chainsSTPtr->qryArySI[slElm] = -1;
         chainsSTPtr->chainArySI[slElm] = 0;
         chainsSTPtr->scoreArySI[slElm] = 0;
         chainsSTPtr->gapScoreArySI[slElm] = 0;
         chainsSTPtr->nextArySI[slElm++] = 0;
      } /*Loop: blank chain arrays*/

   } /*If: need to blank chain arrays*/

   if(chainsSTPtr->indexArySI)
   { /*If: need to blank chain arrays*/

      slElm = 0;

      while(slElm < chainsSTPtr->lenIndexSI)
         chainsSTPtr->indexArySI[slElm++] = -1;

   } /*If: need to blank chain arrays*/

   chainsSTPtr->lenIndexSI = 0;
   chainsSTPtr->lenChainSI = 0;
   chainsSTPtr->lastIndexSI = 0;
} /*blank_chains_mapRead*/

/*-------------------------------------------------------\
| Fun02: init_chains_mapRead
|   - initialize a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct to initialize
| Output:
|   - Modifies:
|     o refKmerAryUI, qryKmerAryUI, chainLenArySI,
|       scoreArySI, and indexArySI to be 0 (null)
|     o lenChainSI, sizeChainSI, lenIndexSL, and
|       sizeIndexSI in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
init_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to initialize*/
){
   if(! chainsSTPtr)
      return; /*nothing to initialize*/

   chainsSTPtr->refArySI = 0;
   chainsSTPtr->qryArySI = 0;
   chainsSTPtr->chainArySI = 0;
   chainsSTPtr->scoreArySI = 0;
   chainsSTPtr->gapScoreArySI = 0;
   chainsSTPtr->nextArySI = 0;
   chainsSTPtr->indexArySI = 0;

   chainsSTPtr->sizeChainSI = 0;
   chainsSTPtr->sizeIndexSI = 0;

   blank_chains_mapRead(chainsSTPtr);
} /*init_chains_mapReads*/

/*-------------------------------------------------------\
| Fun03: freeStack_chains_mapRead
|   - frees all variables in a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead with variables to free
| Output:
|   - Modifies:
|     o frees refKmerAryUI, qryKmerAryUI, chainLenArySI,
|       scoreArySI, and indexArySI and sets to 0 (null)
|     o lenChainSI, sizeChainSI, lenIndexSL, and
|       sizeIndexSI in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
freeStack_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to free vars*/
){
   if(! chainsSTPtr)
      return;

   if(chainsSTPtr->refArySI)
      free(chainsSTPtr->refArySI);

   if(chainsSTPtr->qryArySI)
      free(chainsSTPtr->qryArySI);

   if(chainsSTPtr->chainArySI)
      free(chainsSTPtr->chainArySI);

   if(chainsSTPtr->scoreArySI)
      free(chainsSTPtr->scoreArySI);

   if(chainsSTPtr->gapScoreArySI)
      free(chainsSTPtr->gapScoreArySI);

   if(chainsSTPtr->nextArySI)
      free(chainsSTPtr->nextArySI);

   if(chainsSTPtr->indexArySI)
      free(chainsSTPtr->indexArySI);

   init_chains_mapRead(chainsSTPtr);
} /*freeStack_chains_mapRead*/

/*-------------------------------------------------------\
| Fun04: freeHeap_chains_mapRead
|   - frees a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to free
| Output:
|   - Modifies:
|     o frees chainsSTPtr (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to free*/
){
   if(! chainsSTPtr)
      return;

   freeStack_chains_mapRead(chainsSTPtr);
   free(chainsSTPtr);
} /*freeHeap_chains_mapRead*/

/*-------------------------------------------------------\
| Fun05: setup_chains_mapRead
|   - sets up memory for a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to setup memory
|   - lenChainsSL:
|     o length of chains arrays
|   - lenRefUI:
|     o length of reference
| Output:
|   - Modifies:
|     o allocates lenChainsSL items of memory to
|       refCkmerArySI, qryKmerArySI, chainLenArySI,
|       and scoreArySI in chainsSTPtr
|     o allocates lenRefUI memory to indexArySI in
|       chainsSTPtr
|     o sizeChainSI in chainsSTPtr to be lenChainSI 
|     o sizeIndexSI in chainsSTPtr to be lenIndexSL 
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*to setup*/
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
){
   signed long slElm = 0;

   if(! chainsSTPtr)
      return 0;

   /*making sure no old memory hanging around*/
   freeStack_chains_mapRead(chainsSTPtr);

   chainsSTPtr->refArySI =
      malloc(lenChainsSL * sizeof(signed int));
   if(! chainsSTPtr->refArySI)
      goto memErr_fun05;

   chainsSTPtr->qryArySI =
      malloc(lenChainsSL * sizeof(signed int));
   if(! chainsSTPtr->qryArySI)
      goto memErr_fun05;

   chainsSTPtr->chainArySI =
      calloc(
         lenChainsSL,
         sizeof(signed int)
      );
   if(! chainsSTPtr->chainArySI)
      goto memErr_fun05;

   chainsSTPtr->scoreArySI =
      calloc(
         lenChainsSL,
         sizeof(signed int)
      );
   if(! chainsSTPtr->scoreArySI)
      goto memErr_fun05;

   chainsSTPtr->gapScoreArySI =
      calloc(
         lenChainsSL,
         sizeof(signed int)
      );
   if(! chainsSTPtr->gapScoreArySI)
      goto memErr_fun05;

   chainsSTPtr->nextArySI =
      calloc(
         lenChainsSL,
         sizeof(signed int)
      );
   if(! chainsSTPtr->nextArySI)
      goto memErr_fun05;


   chainsSTPtr->indexArySI =
      calloc(
         lenRefUI,
         sizeof(signed int)
      );
   if(! chainsSTPtr->indexArySI)
      goto memErr_fun05;


   chainsSTPtr->sizeChainSI = lenChainsSL;
   chainsSTPtr->sizeIndexSI = (signed int) lenRefUI;


   slElm = 0;
   while(slElm < lenChainsSL)
   { /*Loop: blank kmer starting coords arrays*/
      chainsSTPtr->refArySI[slElm] = -1;
      chainsSTPtr->qryArySI[slElm++] = -1;
   } /*Loop: blank kmer starting coords arrays*/

   slElm = 0;
   while(slElm < lenRefUI)
      chainsSTPtr->indexArySI[slElm++] = -1;


   return 0;

   memErr_fun05:;
      return def_memErr_mapRead;
} /*setup_chains_mapRead*/

/*-------------------------------------------------------\
| Fun06: mk_chains_mapRead
|   - makes a heap allocated chains_mapRead struct
| Input:
|   - lenChainsSL:
|     o length of chains arrays
|   - lenRefUI:
|     o length of reference
| Output:
|   - Returns:
|     o pointer to initialized and setup chains_mapRead
|       struct
|     o 0 for memory errors
\-------------------------------------------------------*/
struct chains_mapRead *
mk_chains_mapRead(
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
){
   struct chains_mapRead *retSTPtr = 0;

   retSTPtr = malloc(sizeof(struct chains_mapRead));
   if(! retSTPtr)
      goto memErr_fun06;

   init_chains_mapRead(retSTPtr);

   if(
      setup_chains_mapRead(
         retSTPtr,
         lenChainsSL,
         (signed long) lenRefUI
      )
   ) goto memErr_fun06;


   goto ret_fun06;

   memErr_fun06:;
      if(retSTPtr)
         freeHeap_chains_mapRead(retSTPtr);
      retSTPtr = 0;
      goto ret_fun06;
   ret_fun06:;
      return retSTPtr;
} /*mk_chains_mapRead*/

/*-------------------------------------------------------\
| Fun07: realloc_chains_mapRead
|   - reallocates memroy in a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to rellocate memory
|   - lenChainsSL:
|     o new length of chains arrays (0 for do not resize)
|   - lenRefUI:
|     o new length of reference (0 for do not resize)
| Output:
|   - Modifies:
|     o refKmerArySI, qryKmerArySI, chainLenArySI,
|       scoreArySI in chainsSTPtr to be reallocated to
|       lenChainSI if lenChainSI >chainsSTPtr->sizeChainSI
|     o indexArySI chainsSTPtr to be reallocated to
|       lenRefUI if lenRefUI  > chainsSTPtr->sizeIndexSI
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
realloc_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*reallocate*/
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC:
   '   - reallocates memroy in a chains_mapRead struct
   '   o fun07 sec01:
   '     - variable declaration + input check
   '   o fun07 sec02:
   '     - chain arrays reallocs + initialize
   '   o fun07 sec03:
   '     - reference index array rellocs
   '   o fun07 sec04:
   '     - return error or no error
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^   - variable declaration + input check
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *scPtr = 0;
   signed long slElm = 0;

   if(! chainsSTPtr)
      return 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^   - chain arrays reallocs + initialize
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(lenChainsSL > chainsSTPtr->sizeChainSI)
   { /*If: reallocating chains array memory*/
      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->refArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->refArySI = (signed int *) scPtr;

      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->qryArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->qryArySI = (signed int *) scPtr;

      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->chainArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->chainArySI = (signed int *) scPtr;

      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->scoreArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->scoreArySI = (signed int *) scPtr;

      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->gapScoreArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->gapScoreArySI = (signed int *) scPtr;


      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->nextArySI,
            lenChainsSL * sizeof(signed int)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->nextArySI = (signed int *) scPtr;


      slElm = chainsSTPtr->sizeChainSI;
      while(slElm < lenChainsSL)
      { /*Loop: blank kmer starting coords arrays*/
         chainsSTPtr->refArySI[slElm] = -1;
         chainsSTPtr->qryArySI[slElm] = -1;
         chainsSTPtr->chainArySI[slElm] = 0;
         chainsSTPtr->scoreArySI[slElm] = 0;
         chainsSTPtr->gapScoreArySI[slElm] = 0;
         chainsSTPtr->nextArySI[slElm++] = 0;
      } /*Loop: blank kmer starting coords arrays*/

      chainsSTPtr->sizeChainSI = lenChainsSL;
   } /*If: reallocating chains array memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^   - reference index array rellocs
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(lenRefUI > (unsigned int) chainsSTPtr->sizeIndexSI)
   { /*If: reallocating chains array memory*/
      scPtr =
         (signed char *)
         realloc(
            chainsSTPtr->indexArySI,
            lenRefUI * sizeof(signed long)
         );
      if(! scPtr)
         goto memErr_fun07_sec04;
      chainsSTPtr->indexArySI = (signed int *) scPtr;

      slElm = chainsSTPtr->sizeIndexSI;
      while(slElm < lenRefUI)
         chainsSTPtr->indexArySI[slElm++] = 0;

      chainsSTPtr->sizeIndexSI = lenRefUI;
   } /*If: reallocating chains array memory*/


   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec04:
   ^   - return error or no error
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun07_sec04:;
      return def_memErr_mapRead;
} /*realloc_chains_mapRead*/

/*-------------------------------------------------------\
| Fun08: blank_set_mapRead
|   - blank a set_mapRead struct
| Input:
|   - setSTPtr:
|     o pointer to set_mapRead struct to blank
| Output:
|   - Modifies:
|     o sets values in setSTPtr to defaults
\-------------------------------------------------------*/
void
blank_set_mapRead(
   struct set_mapRead *setSTPtr
){
   if(! setSTPtr)
      return;

   setSTPtr->alnEndsBl = def_alnEnd_defsMapRead;
   setSTPtr->subBl = def_subAln_defsMapRead;

   setSTPtr->gapSI = def_chainGap_defsMapRead;
   setSTPtr->matchSI = def_chainMatch_defsMapRead;

   setSTPtr->minScoreF = def_minScore_defsMapRead;
   setSTPtr->minMatchF = def_maxPercLen_defsMapRead;
   setSTPtr->minPercLenF = def_minPercLen_defsMapRead;
   setSTPtr->chainMinNtF = def_chainMinLen_defsMapRead;
   setSTPtr->maxLenPercF = def_maxPercLen_defsMapRead;

   setSTPtr->lenArySI[0] = def_minLen_defsMapRead;
   setSTPtr->lenArySI[1] = def_secLen_defsMapRead;
   setSTPtr->lenArySI[2] = def_thirdLen_defsMapRead;

   setSTPtr->kmerAryUC[0] = def_minKmer_defsMapRead;
   setSTPtr->kmerAryUC[1] = def_secKmer_defsMapRead;
   setSTPtr->kmerAryUC[2] = def_thirdKmer_defsMapRead;

   setSTPtr->chainLenArySI[0] = def_minChain_defsMapRead;
   setSTPtr->chainLenArySI[1] = def_secChain_defsMapRead;
   setSTPtr->chainLenArySI[2] =
        def_thirdChain_defsMapRead;

   setSTPtr->lenKmersUC = 3;

   if(setSTPtr->alnSetST)
      init_alnSet(setSTPtr->alnSetST);
      /*alnSet has no blank function*/
} /*blank_set_mapRead*/

/*-------------------------------------------------------\
| Fun09: init_set_mapRead
|   - initializes a set_mapRead struct
| Input:
|   - setSTPtr:
|     o pointer to set_mapRead struct to initializes
| Output:
|   - Modifies:
|     o sets variables to defaults and pointers to null
\-------------------------------------------------------*/
void
init_set_mapRead(
   struct set_mapRead *setSTPtr
){
   if(! setSTPtr)
      return;

   setSTPtr->alnSetST = 0;

   blank_set_mapRead(setSTPtr);
} /*init_set_mapRead*/

/*-------------------------------------------------------\
| Fun10: freeStack_set_mapRead
|   - frees variables in a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer with variables to free
| Output:
|   - Modifies:
|     o frees all pointers and intitializes; defaults/null
\-------------------------------------------------------*/
void
freeStack_set_mapRead(
   struct set_mapRead *setSTPtr
){
   if(! setSTPtr)
      return;

   if(setSTPtr->alnSetST)
      freeHeap_alnSet(setSTPtr->alnSetST);

   init_set_mapRead(setSTPtr);
} /*freeStack_set_mapRead*/

/*-------------------------------------------------------\
| Fun11: freeHeap_set_mapRead
|   - frees a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o frees setSTPtr (you must set to null/0)
\-------------------------------------------------------*/
void
freeHeap_set_mapRead(
   struct set_mapRead *setSTPtr
){
   if(! setSTPtr)
      return;

   freeStack_set_mapRead(setSTPtr);
   free(setSTPtr);
} /*freeHeap_set_mapRead*/

/*-------------------------------------------------------\
| Fun12: setup_set_mapRead
|   - allocates memory for structs in a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o alnSetST to be on heap
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_set_mapRead(
   struct set_mapRead *setSTPtr
){
   if(! setSTPtr)
      return 0;

   if(! setSTPtr->alnSetST)
   { /*If: do not have alignment settings*/

      setSTPtr->alnSetST = malloc(sizeof(struct alnSet));
      if(! setSTPtr->alnSetST)
         goto memErr_fun12;
      init_alnSet(setSTPtr->alnSetST);

   } /*If: do not have alignment settings*/

   return 0;

   memErr_fun12:;
      return def_memErr_mapRead;
} /*setup_set_mapRead*/

/*-------------------------------------------------------\
| Fun13: blank_ref_mapRead
|   - blanks (sets to defaults) a ref_mapRead structure
|   - use addRef_ref_mapRead to change references
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to blank
| Output:
|   - Modifies:
|     o sets all values in refSTPtr to defaults or for
|       arrays null
\-------------------------------------------------------*/
void
blank_ref_mapRead(
   struct ref_mapRead *refSTPtr
){
   signed int iterSI = 0;

   if(! refSTPtr)
      return;

   /*these are not good defaults, but then they also will
   `  be set when the user adds the reference sequence
   */
   refSTPtr->lenKmerUC = def_thirdKmer_defsMapRead;
   refSTPtr->minChainLenSI = def_thirdChain_defsMapRead;

   /*need to blank all kmers in arrays*/
   for(
      iterSI = 0;
      iterSI < refSTPtr->lenSI;
      ++iterSI
   ){ /*Loop: blank kmer arrays*/
      if(refSTPtr->kmerArySI)
        refSTPtr->kmerArySI[0] = def_endKmers_kmerFun;
      if(refSTPtr->indexArySI)
        refSTPtr->indexArySI[0] = (signed int) iterSI;
   }  /*Loop: blank kmer arrays*/

   refSTPtr->lenSI = 0;

   if(refSTPtr->seqSTPtr)
      blank_seqST(refSTPtr->seqSTPtr);
} /*blank_ref_mapRead*/

/*-------------------------------------------------------\
| Fun14: init_ref_mapRead
|   - initialize a ref_mapRead structure
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to initialize
| Output:
|   - Modifies:
|     o sets pointers to 0/null & other values to default
\-------------------------------------------------------*/
void
init_ref_mapRead(
   struct ref_mapRead *refSTPtr
){
   if(! refSTPtr)
     return;

   refSTPtr->seqSTPtr = 0;
   refSTPtr->kmerArySI = 0;
   refSTPtr->indexArySI = 0;
   refSTPtr->sizeKmerSI = 0;
   refSTPtr->sizeIndexSI = 0;
   refSTPtr->lenSI = 0;

   blank_ref_mapRead(refSTPtr);
} /*init_ref_mapRead*/

/*-------------------------------------------------------\
| Fun15: freeStack_ref_mapRead
|   - frees variables in a ref_mapRead struct
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer with variables to free
| Output:
|   - Modifies:
|     o frees all heap allocated variabes in refSTPtr
|     o all variabes in refSTPtr to be 0/null (free) or
|       defaults
\-------------------------------------------------------*/
void
freeStack_ref_mapRead(
   struct ref_mapRead *refSTPtr
){
   if(! refSTPtr)
      return;

    if(refSTPtr->seqSTPtr)
       freeHeap_seqST(refSTPtr->seqSTPtr);
    refSTPtr->seqSTPtr = 0;

    if(refSTPtr->kmerArySI)
       free(refSTPtr->kmerArySI);
    if(refSTPtr->indexArySI)
       free(refSTPtr->indexArySI);

    init_ref_mapRead(refSTPtr);
} /*freeStack_ref_mapRead*/

/*-------------------------------------------------------\
| Fun16: freeHeap_ref_mapRead
|   - frees a ref_mapRead struct
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o frees refSTPtr (you must set to 0/null)
\-------------------------------------------------------*/
void
freeHeap_ref_mapRead(
   struct ref_mapRead *refSTPtr
){
   if(! refSTPtr)
      return;
   freeStack_ref_mapRead(refSTPtr);
   free(refSTPtr);
} /*freeHeap_ref_mapRead*/

/*-------------------------------------------------------\
| Fun17: addRef_ref_mapRead
|   - adds a reference sequence to a ref_map read struct
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to add reference to
|   - seqSTPtr:
|     o seqST struct pointer with reference sequence
|     o must be heap allocated; otherwise messes up
|       freeStack_ref_mapRead
|     o on success, you should set this to null, but do
|       not free
|     o on failure, free seqSTPtr; only assigned on succes
|   - setSTPtr:
|     o set_mapRead struct with settings for alignment
| Output:
|   - Modifies:
|     o seqSTPtr in refSTPtr to point to seqSTPtr
|     o on success sequence in seqSTPtr is converted
|       to alnSet lookup index's
|     o kmerArySI to have refernce kmers
|     o sortKmerArySI to have refernce kmers and be sorted
|     o indexArySI to have index's to convert
|       sortKmerArySI index to kmerArySI index
|     o lenKmerUC and minChainLenSI to bet set to their
|       correct values in setSTPtr
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|     o def_noStruct_mapRead if any input is 0/null
\-------------------------------------------------------*/
signed char
addRef_ref_mapRead(
   struct ref_mapRead *refSTPtr, /*add reference to*/
   struct seqST *seqSTPtr,       /*reference to add*/
   struct set_mapRead *setSTPtr  /*has settings to add*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun17 TOC:
   '   - adds a reference sequence to ref_map read struct
   '   o fun17 sec01:
   '     - variable declaration and make sure have input
   '   o fun17 sec02:
   '     - find reference alignment settings
   '   o fun17 sec03:
   '     - build kmer arrays
   '   o fun17 sec04:
   '     - allocate memory for sorted and index arrays
   '   o fun17 sec05:
   '     - copy kmers into sorted array and sort kmers
   '   o fun17 sec06:
   '     - set seqsTPtr pointer and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec01:
   ^   - variable declaration and make sure have input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int siKmer = 0;
   signed int siPos = 0;

   if(! refSTPtr)
      goto noStruct_fun17_sec06;
   else if(! seqSTPtr)
      goto noStruct_fun17_sec06;
   else if(! seqSTPtr->seqStr)
      goto noStruct_fun17_sec06;
   else if(seqSTPtr->seqStr[0] == '\0')
      goto noStruct_fun17_sec06;
   else if(! setSTPtr)
      goto noStruct_fun17_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec02:
   ^   - find reference alignment settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(
      siKmer = 0;
      siKmer < (signed int) setSTPtr->lenKmersUC;
      ++siKmer
   ){ /*Loop: find kmer and chain size*/
      if(
            refSTPtr->seqSTPtr->seqLenSL
          < setSTPtr->lenArySI[siKmer]
      ) break; /*found kmer settings*/
   }  /*Loop: find kmer and chain size*/

   siKmer -=
      (siKmer == (signed int) setSTPtr->lenKmersUC);
      /*subtracts if hit end of kmer array (max length)*/
  
   refSTPtr->lenKmerUC = setSTPtr->kmerAryUC[siKmer];
   refSTPtr->minChainLenSI =
      setSTPtr->chainLenArySI[siKmer];

   /*needed indexed bases for kmer conversion*/
   seqToIndex_alnSet(seqSTPtr->seqStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec03:
   ^   - build kmer arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   refSTPtr->lenSI =
      seqToKmer_kmerFun(
         seqSTPtr->seqStr,
         seqSTPtr->seqLenSL,
         &refSTPtr->kmerArySI,
         &refSTPtr->sizeKmerSI,
         refSTPtr->lenKmerUC
      );

   if(! refSTPtr->lenSI)
      goto memErr_fun17_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec04:
   ^   - allocate memory for sorted and index arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(refSTPtr->sizeIndexSI < refSTPtr->lenSI)
   { /*If: need more memory*/

      if(refSTPtr->indexArySI)
         free(refSTPtr->indexArySI);
      refSTPtr->indexArySI = 0;

      refSTPtr->indexArySI =
         malloc(
            (refSTPtr->lenSI + 8) * sizeof(signed int)
         );

      if(! refSTPtr->indexArySI)
         goto memErr_fun17_sec06;

      refSTPtr->sizeIndexSI = refSTPtr->lenSI;
   } /*If: need more memory*/

   for(siPos = 0; siPos < refSTPtr->lenSI; ++siPos)
      refSTPtr->indexArySI[siPos] = siPos;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec05:
   ^   - copy kmers into sorted array and sort kmers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uiTwinSort_shellSort(
      (unsigned int *) refSTPtr->kmerArySI,
      (unsigned int *) refSTPtr->indexArySI,
      0,
      refSTPtr->lenSI - 1
   ); /*using unsigned int to force -1's to end of array*/

   refSTPtr->lenSI =
      seqToKmer_kmerFun(
         seqSTPtr->seqStr,
         seqSTPtr->seqLenSL,
         &refSTPtr->kmerArySI,
         &refSTPtr->sizeKmerSI,
         refSTPtr->lenKmerUC
      );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun17 Sec06:
   ^   - set seqSTPtr pointer and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   refSTPtr->seqSTPtr = seqSTPtr;
   return 0;

   memErr_fun17_sec06:;
     siKmer = def_memErr_mapRead;
     goto errRet_fun17_sec06;

   noStruct_fun17_sec06:;
     siKmer = def_noStruct_mapRead;
     goto errRet_fun17_sec06;

   errRet_fun17_sec06:;
      indexToSeq_alnSet(seqSTPtr->seqStr);
      return (signed char) siKmer;
} /*addRef_ref_mapRead*/

/*-------------------------------------------------------\
| Fun18: mk_ref_mapRead
|   - makes and initializes a ref_mapRead structure
| Input:
|   - seqSTPtr:
|     o seqST struct pointer with reference sequence
|     o must be heap allocated; otherwise messes up
|       freeStack_ref_mapRead
|     o on success, you should set this to null, but do
|       not free
|     o on failure, free seqSTPtr; only assigned on succes
|   - setSTPtr:
|     o set_mapRead struct with settings for alignment
| Output:
|   - Returns:
|     o 0 for errors
|     o pointer to setup ref_mapRead structure
\-------------------------------------------------------*/
struct ref_mapRead *
mk_ref_mapRead(
   struct seqST *seqSTPtr,       /*reference to add*/
   struct set_mapRead *setSTPtr  /*has settings to add*/
){
   struct ref_mapRead *refSTPtr = 0;
  
   if(! seqSTPtr)
      goto err_fun18;
   else if(! seqSTPtr->seqStr)
      goto err_fun18;
   else if(! seqSTPtr->seqStr[0])
      goto err_fun18;
   else if(! setSTPtr)
      goto err_fun18;

   refSTPtr = malloc( sizeof(struct ref_mapRead) );
   if(! refSTPtr)
      goto err_fun18;

   init_ref_mapRead(refSTPtr);

   if(
      addRef_ref_mapRead(
         refSTPtr,
         seqSTPtr,
         setSTPtr
      )
   ) goto err_fun18;

   goto ret_fun18;

   err_fun18:;
      if(refSTPtr)
         freeHeap_ref_mapRead(refSTPtr);
      refSTPtr = 0;
      goto ret_fun18;

   ret_fun18:;
      return refSTPtr;
} /*mk_ref_mapRead*/

/*-------------------------------------------------------\
| Fun19: getRef_ref_mapRead
|   - gets and stores next reference in ref_mapRead struct
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to store reference in
|   - setSTPtr:
|     o set_mapRead struct with settings for alignment
|   - fxBl:
|     o tells if file is fastq (1) or regular fasta (0)
|   - gzBl:
|     o tells if file is gziped (1) or regular fasta (0)
|     o TODO: add gzip support
|   - refFILE:
|     o FILE pointer to file with reference to store
| Output:
|   - Modifies:
|     o seqSTPtr in refSTPtr to have reference sequence
|     o kmerArySI in refSTPtr to have refernce kmers
|     o sortKmerArySI to have refernce kmers and be sorted
|     o indexArySI to have index's to convert
|       sortKmerArySI index to kmerArySI index
|     o lenKmerUC and minChainLenSI to bet set to their
|       correct values in setSTPtr
|     o refFILE to be on next entry in fastx file
|   - Returns:
|     o 0 for errors
|     o def_EOF_mapRead if hit end of file
|     o def_memErr_mapRead for memory errors
|     o def_noStruct_mapRead if setSTPtr or refSTPtr is 0
|     o def_badFile_mapRead for invalid fastx entries
|       * also for gz files (currently not supported)
\-------------------------------------------------------*/
signed char
getRef_ref_mapRead(
   struct ref_mapRead *refSTPtr,/*gets reference seq*/
   struct set_mapRead *setSTPtr,/*has settings to add*/
   signed char fxBl,            /*1: fastq; 0: fasta*/
   signed char gzBl,            /*1: gzipped; 0: not*/
   void *refFILE                /*fx file with reference*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun19 TOC:
   '   - gets and stores next reference in ref_mapRead
   '   o fun19 sec01:
   '     - variable declarations
   '   o fun19 sec02:
   '     - checks and memory allocation
   '   o fun19 sec03:
   '     - get and add reference sequence
   '   o fun19 sec04:
   '     - setup reference in refSTPtr
   '   o fun19 sec05:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec01:
   ^   - variable declarations
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   signed char errSC = 0;
   signed char eofSC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec02:
   ^   - checks and memory allocation
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   if(! refSTPtr)
      goto noStruct_fun19_sec05;
   else if(! setSTPtr)
      goto noStruct_fun19_sec05;
   else if(gzBl)
      goto badFile_fun19_sec05;
      /*currently not supported*/

   if(! refSTPtr->seqSTPtr)
   { /*If: need a sequence structure*/
      refSTPtr->seqSTPtr = malloc( sizeof(struct seqST) );

      if(! refSTPtr->seqSTPtr)
         goto memErr_fun19_sec05;
      init_seqST(refSTPtr->seqSTPtr);
   } /*If: need a sequence structure*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec03:
   ^   - get and add reference sequence
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   /*TODO add gzip support*/
   if(fxBl)
   { /*If: fastq file*/
      errSC = getFq_seqST(refFILE, refSTPtr->seqSTPtr);
   } /*If: fastq file*/

   else
   { /*Else: fasta file*/
      errSC = getFa_seqST(refFILE, refSTPtr->seqSTPtr);
   } /*Else: fasta file*/

   if(errSC)
   { /*If: had error*/
      if(errSC & def_fileErr_seqST)
         goto badFile_fun19_sec05;
      else if(errSC & def_memErr_seqST)
         goto memErr_fun19_sec05;
      else if(refSTPtr->seqSTPtr->seqStr[0] == '\0')
         goto eof_fun19_sec05;
         /*end of file*/
      else
         eofSC = 1; /*last sequence in file*/
   } /*If: had error*/
  
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec04:
   ^   - setup reference in refSTPtr
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   errSC =
      addRef_ref_mapRead(
         refSTPtr,
         refSTPtr->seqSTPtr,
         setSTPtr
      );

   if(errSC)
      goto memErr_fun19_sec05;
      /*already know structures exist, so must be memory
      `   error
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec05:
   ^   - return result
   \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

   if(eofSC)
      errSC = def_EOF_mapRead;
   else
      errSC = 0;

   goto ret_fun19_sec05;

   memErr_fun19_sec05:;
      errSC = def_memErr_mapRead;
      goto ret_fun19_sec05;

   noStruct_fun19_sec05:;
      errSC = def_noStruct_mapRead;
      goto ret_fun19_sec05;

   badFile_fun19_sec05:;
      errSC = def_badFile_mapRead;
      goto ret_fun19_sec05;

   eof_fun19_sec05:;
      errSC = def_EOF_mapRead;
      goto ret_fun19_sec05;

   ret_fun19_sec05:;
      return errSC;
} /*getRef_ref_mapRead*/

/*-------------------------------------------------------\
| Fun20: getRefBin_ref_mapRead
|   - gets a binary file with reference
| Input:
|   - refSTPtr:
|     o ref_mapRead struct to add reference to
|   - refFILE:
|     o FILE pointer to binary file with reference
|       sequence, kmers, sorted index's, minimum chain
|       length, and kmer length
|     o Format:
|       1: kmer length [of one kmer] (unsigned char)
|       2: minimum kmers in a chain (signed int)
|       3: length of reference id (name) (signed int)
|       4: reference id (signed char)
|       5: length of reference sequence (signed int)
|       6: reference sequence (signed char)
|       7: length of kmer sequence (signed int)
|       8: reference kmer sequence (signed int)
|       9: sorted index's for kmer sequence (signed int)
| Output:
|   - Modifies:
|     o refSTPtr to have values in refFILE
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|     o def_noStruct_mapRead if no structure input
|     o def_badFile_mapRead if no file input or if had
|       file error
\-------------------------------------------------------*/
signed char
getRefBin_ref_mapRead(
   struct ref_mapRead *refSTPtr,/*gets reference*/
   void *refFILE               /*bin file with reference*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun20 TOC:
   '   - gets a binary file with reference
   '   o fun20 sec01:
   '     - variable declartions
   '   o fun20 sec02:
   '     - check input & allocate memory for seqST struct
   '   o fun20 sec03:
   '     - get kmer size and minimum chain length
   '   o fun20 sec04:
   '     - get reference id
   '   o fun20 sec05:
   '     - get reference sequence
   '   o fun20 sec06:
   '     - get kmer sequence and sorted kmer index's
   '   o fun20 sec07:
   '     - return result/errors
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec01:
   ^   - variable declartions
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned long tmpUL = 0;
   signed int tmpSI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec02:
   ^   - check input and allocate memory for seqST struct
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! refSTPtr)
      goto noStruct_fun20_sec07;
   if(! refFILE)
      goto noFile_fun20_sec07;

   if(! refSTPtr->seqSTPtr)
   { /*If: need new sequence struct*/
      refSTPtr->seqSTPtr = malloc(sizeof(struct seqST));
      if(! refSTPtr->seqSTPtr)
         goto memErr_fun20_sec07;
      init_seqST(refSTPtr->seqSTPtr);
   } /*If: need new sequence struct*/

   else
      blank_seqST(refSTPtr->seqSTPtr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec03:
   ^   - get kmer size and minimum chain length
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tmpUL =
      fread(
         &refSTPtr->lenKmerUC,
         sizeof(unsigned char),
         1,
         (FILE *) refFILE
      );

   if(! tmpUL)
      goto badEntry_fun20_sec07;

   tmpUL =
      fread(
         &refSTPtr->minChainLenSI,
         sizeof(signed int),
         1,
         (FILE *) refFILE
      );

   if(! tmpUL)
      goto badEntry_fun20_sec07;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec04:
   ^   - get reference id
   ^   o fun20 sec04 sub01:
   ^     - get id length from file
   ^   o fun20 sec04 sub02:
   ^     - allocate memory for id string
   ^   o fun20 sec04 sub03:
   ^     - get id from file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun20 Sec04 Sub01:
   *   - get id length from file
   \*****************************************************/

   tmpUL =
      fread(
         &refSTPtr->lenSI,
         sizeof(signed int),
         1,
         (FILE *) refFILE
      );

   if(! tmpUL)
      goto badEntry_fun20_sec07;

   /*****************************************************\
   * Fun20 Sec04 Sub02:
   *   - allocate memory for id string
   \*****************************************************/

   if(refSTPtr->seqSTPtr->idStr)
   { /*If: have id buffer*/

      if(refSTPtr->seqSTPtr->idSizeSL < refSTPtr->lenSI)
      { /*If: id buffer is to short*/
         free(refSTPtr->seqSTPtr->idStr);
         refSTPtr->seqSTPtr->idStr = 0;
      }  /*If: id buffer is to short*/
   } /*If: have id buffer*/

   if(! refSTPtr->seqSTPtr->idStr)
   { /*If: need more memory*/
      refSTPtr->seqSTPtr->idStr =
         malloc(
           (refSTPtr->lenSI + 8) * sizeof(signed char)
         );
      if(! refSTPtr->seqSTPtr->idStr)
         goto memErr_fun20_sec07;
      refSTPtr->seqSTPtr->idSizeSL = refSTPtr->lenSI;
   } /*If: need more memory*/

   /*****************************************************\
   * Fun20 Sec04 Sub03:
   *   - get id from file
   \*****************************************************/

   refSTPtr->seqSTPtr->idLenSL =
      fread(
         refSTPtr->seqSTPtr->idStr,
         sizeof(signed char),
         refSTPtr->lenSI,
         (FILE *) refFILE
      );

   if(! refSTPtr->seqSTPtr->idLenSL)
      goto badEntry_fun20_sec07;
   if(refSTPtr->seqSTPtr->idLenSL < refSTPtr->lenSI)
      goto badEntry_fun20_sec07; /*hit EOF early*/

   refSTPtr->seqSTPtr->idStr[
      refSTPtr->seqSTPtr->idLenSL
   ] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec05:
   ^   - get reference sequence
   ^   o fun20 sec05 sub01:
   ^     - get sequence length from file
   ^   o fun20 sec05 sub02:
   ^     - allocate memory for sequence
   ^   o fun20 sec05 sub03:
   ^     - set q-score entry to null
   ^   o fun20 sec05 sub04:
   ^     - get sequence from file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun20 Sec05 Sub01:
   *   - get sequence length from file
   \*****************************************************/

   tmpUL =
      fread(
         &refSTPtr->lenSI,
         sizeof(signed int),
         1,
         (FILE *) refFILE
      );

   if(! tmpUL)
      goto badEntry_fun20_sec07;

   /*****************************************************\
   * Fun20 Sec05 Sub02:
   *   - allocate memory for sequence
   \*****************************************************/

   if(refSTPtr->seqSTPtr->seqStr)
   { /*If: have sequence buffer*/

      if(refSTPtr->seqSTPtr->seqSizeSL < refSTPtr->lenSI)
      { /*If: need more memory*/
         free(refSTPtr->seqSTPtr->seqStr);
         refSTPtr->seqSTPtr->seqStr = 0;
      }  /*If: need more memory*/
   } /*If: have sequence buffer*/

   if(! refSTPtr->seqSTPtr->seqStr)
   { /*If: need to allocate memory*/
      refSTPtr->seqSTPtr->seqStr =
         malloc(
            (refSTPtr->lenSI + 8) * sizeof(signed char)
         );
      if(! refSTPtr->seqSTPtr->seqStr)
         goto memErr_fun20_sec07;
      refSTPtr->seqSTPtr->seqSizeSL = refSTPtr->lenSI;
   } /*If: need to allocate memory*/

   /*****************************************************\
   * Fun20 Sec05 Sub03:
   *   - set q-score entry to null
   \*****************************************************/

   if(! refSTPtr->seqSTPtr->qStr)
   { /*If: need memory for q-score entry*/
      refSTPtr->seqSTPtr->qStr =
         malloc(9 * sizeof(signed char));
      if(! refSTPtr->seqSTPtr->qStr)
         goto memErr_fun20_sec07;
      refSTPtr->seqSTPtr->qSizeSL = 1;
   } /*If: need memory for q-score entry*/

   refSTPtr->seqSTPtr->qStr[0] = '\0';
   refSTPtr->seqSTPtr->qLenSL = 0;

   /*****************************************************\
   * Fun20 Sec05 Sub04:
   *   - get sequence from file
   \*****************************************************/

   refSTPtr->seqSTPtr->seqLenSL =
      fread(
         refSTPtr->seqSTPtr->seqStr,
         sizeof(signed char),
         refSTPtr->lenSI,
         (FILE *) refFILE
      );

   if(! refSTPtr->seqSTPtr->seqLenSL)
      goto badEntry_fun20_sec07;
   if(
        refSTPtr->seqSTPtr->seqLenSL
      < (unsigned int) refSTPtr->lenSI
   ) goto badEntry_fun20_sec07;
     /*hit EOF early*/

   refSTPtr->seqSTPtr->seqStr[
      refSTPtr->seqSTPtr->seqLenSL
   ] = '\0';

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec06:
   ^   - get kmer sequence and sorted kmer index's
   ^   o fun20 sec06 sub01:
   ^     - get kmer sequence length
   ^   o fun20 sec06 sub03:
   ^     - allocate memory for kmer sequence
   ^   o fun20 sec06 sub04:
   ^     - allocate memory for sorted kmer indexs
   ^   o fun20 sec06 sub04:
   ^     - get kmer sequence and sorted kmer index
   ^   o fun20 sec06 sub05:
   ^     - get sorted kmer index's
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun20 Sec06:
   *   - get kmer sequence and sorted kmer index's
   \*****************************************************/

   tmpUL =
      fread(
         &tmpSI,
         sizeof(signed int),
         1,
         (FILE *) refFILE
      );

   if(! tmpUL)
      goto badEntry_fun20_sec07;

   /*****************************************************\
   * Fun20 Sec06 Sub02:
   *   - allocate memory for kmer sequence
   \*****************************************************/

   if(refSTPtr->kmerArySI)
   { /*If: have kmer sequence buffer*/

      if(refSTPtr->sizeKmerSI < tmpSI)
      { /*If: need more memory*/
         free(refSTPtr->kmerArySI);
         refSTPtr->kmerArySI = 0;
      }  /*If: need more memory*/
   } /*If: have kmer sequence buffer*/


   if(! refSTPtr->kmerArySI)
   { /*If: need to allocate memory*/
      refSTPtr->kmerArySI =
         malloc((tmpSI + 8) * sizeof(signed int));
      if(! refSTPtr->kmerArySI)
         goto memErr_fun20_sec07;
      refSTPtr->sizeKmerSI = tmpSI;
   } /*If: need to allocate memory*/

   /*****************************************************\
   * Fun20 Sec06 Sub03:
   *   - allocate memory for sorted kmer index's
   \*****************************************************/

   if(refSTPtr->indexArySI)
   { /*If: have kmer sequence buffer*/

      if(refSTPtr->sizeIndexSI < tmpSI)
      { /*If: need more memory*/
         free(refSTPtr->indexArySI);
         refSTPtr->indexArySI = 0;
      }  /*If: need more memory*/
   } /*If: have kmer sequence buffer*/


   if(! refSTPtr->indexArySI)
   { /*If: need to allocate memory*/
      refSTPtr->indexArySI =
         malloc((tmpSI + 8) * sizeof(signed int));
      if(! refSTPtr->indexArySI)
         goto memErr_fun20_sec07;
      refSTPtr->sizeIndexSI = tmpSI;
   } /*If: need to allocate memory*/

   /*****************************************************\
   * Fun20 Sec06 Sub04:
   *   - get kmer sequence
   \*****************************************************/

   tmpUL =
      fread(
         refSTPtr->kmerArySI,
         sizeof(signed int),
         tmpSI,
         (FILE *) refFILE
      );

   if(! tmpUL)
      goto badEntry_fun20_sec07;
   if(tmpUL < (unsigned int) tmpSI)
      goto badEntry_fun20_sec07;
      /*hit EOF early*/

   /*****************************************************\
   * Fun20 Sec06 Sub05:
   *   - get sorted kmer index's
   \*****************************************************/

   refSTPtr->lenSI  =
      fread(
         refSTPtr->indexArySI,
         sizeof(signed int),
         tmpSI,
         (FILE *) refFILE
      );

   if(! refSTPtr->lenSI)
      goto badEntry_fun20_sec07;
   if(refSTPtr->lenSI < tmpSI)
      goto badEntry_fun20_sec07;
      /*hit EOF early*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec07:
   ^   - return result/errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun20_sec07:;
     return def_memErr_mapRead;

   noFile_fun20_sec07:;
     return def_badFile_mapRead;

   noStruct_fun20_sec07:;
     return def_noStruct_mapRead;

   badEntry_fun20_sec07:;
     return def_badFile_mapRead;
} /*getRefBin_ref_mapRead*/

/*-------------------------------------------------------\
| Fun21: writeRefBin_ref_mapRead
|   - writes a ref_mapRead struct to a binary file
| Input:
|   - refSTPtr:
|     o ref_mapRead struct to write
|   - refFILE:
|     o FILE pointer to write contents of refSTPtr to
| Output:
|   - Prints:
|     o refSTPtr to file
|       1: kmer length [of one kmer] (unsigned char)
|       2: minimum kmers in a chain (signed int)
|       3: length of reference id (name) (signed int)
|       4: reference id (signed char)
|       5: length of reference sequence (signed int)
|       6: reference sequence (signed char)
|       7: length of kmer sequence (signed int)
|       8: reference kmer sequence (signed int)
|       9: sorted index's for kmer sequence (signed int)
|   - Returns:
|     o 0 for no errors
|     o def_noStruct_mapRead if refSTPtr is 0/null
|     o def_badFile_mapRead if refFILE is 0/null
\-------------------------------------------------------*/
signed char
writeRefBin_ref_mapRead(
   struct ref_mapRead *refSTPtr,/*reference to write*/
   void *refFILE              /*file to save refSTPtr to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun21 TOC:
   '   - writes a ref_mapRead struct to a binary file
   '   o fun21 sec01:
   '     - variable declarations and initial checks
   '   o fun21 sec02:
   '     - print kmer length and min chain length
   '   o fun21 sec03:
   '     - print id and sequence
   '   o fun21 sec04:
   '     - print kmer sequence and sorted kmer index's
   '   o fun21 sec05:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec01:
   ^   - variable declarations and initial checks
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int tmpSI = 0;

   if(! refSTPtr)
      goto noStruct_fun21_sec05;
   if(! refFILE)
      goto badFile_fun21_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec02:
   ^   - print kmer length and min chain length
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fwrite(
      &refSTPtr->lenKmerUC,
      sizeof(unsigned char),
      1,
      (FILE *) refFILE
   );

   fwrite(
      &refSTPtr->minChainLenSI,
      sizeof(signed int),
      1,
      (FILE *) refFILE
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec03:
   ^   - print id and sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   tmpSI = refSTPtr->seqSTPtr->idLenSL;

   fwrite(
      &tmpSI,
      sizeof(signed int),
      1,
      (FILE *) refFILE
   );

   fwrite(
      refSTPtr->seqSTPtr->idStr,
      sizeof(signed char),
      refSTPtr->seqSTPtr->idLenSL,
      (FILE *) refFILE
   );


   tmpSI = refSTPtr->seqSTPtr->seqLenSL;

   fwrite(
      &tmpSI,
      sizeof(signed int),
      1,
      (FILE *) refFILE
   );

   fwrite(
      refSTPtr->seqSTPtr->seqStr,
      sizeof(signed char),
      refSTPtr->seqSTPtr->seqLenSL,
      (FILE *) refFILE
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec04:
   ^   - print kmer sequence and sorted kmer index's
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fwrite(
      &refSTPtr->lenSI,
      sizeof(signed int),
      1,
      (FILE *) refFILE
   );

   fwrite(
      refSTPtr->kmerArySI,
      sizeof(signed int),
      refSTPtr->lenSI,
      (FILE *) refFILE
   );

   fwrite(
      refSTPtr->indexArySI,
      sizeof(signed int),
      refSTPtr->lenSI,
      (FILE *) refFILE
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec05:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   noStruct_fun21_sec05:;
      return def_noStruct_mapRead;

   badFile_fun21_sec05:;
      return def_badFile_mapRead;
} /*writeRefBin_ref_mapRead*/

/*-------------------------------------------------------\
| Fun22: blank_aln_mapRead
|   - blanks an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer to blank
| Output:
|   - Modifies:
|     o first element kmerArySI in alnSTPtr to be
|       def_endKmers_kmerFun 
|     o lenSI in alnSTPtr to be 0
|     o dirMatrix in alnSTPtr with blank_dirMatrix
|     o chainsSTPtr in alnSTPtr with blank_chains_mapRead
\-------------------------------------------------------*/
void
blank_aln_mapRead(
   struct aln_mapRead *alnSTPtr
){
   if(! alnSTPtr)
      return;

   if(alnSTPtr->kmerArySI)
      alnSTPtr->kmerArySI[0] = def_endKmers_kmerFun;
   alnSTPtr->lenSI = 0;

   if(alnSTPtr->matrixSTPtr)
      blank_dirMatrix(alnSTPtr->matrixSTPtr);

   if(alnSTPtr->chainsSTPtr)
      blank_chains_mapRead(alnSTPtr->chainsSTPtr);
} /*blank_aln_mapRead*/

/*-------------------------------------------------------\
| Fun23: init_aln_mapRead
|   - initializes an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer to initialize
| Output:
|   - Modifies:
|     o kmerArySI, matrixSTPtr, and chainsSTPtr in
|       alnSTPtr to be 0/null
|     o lenSI and sizeSI in alnSTPtr to be 0
\-------------------------------------------------------*/
void
init_aln_mapRead(
   struct aln_mapRead *alnSTPtr
){
   if(! alnSTPtr)
      return;

   alnSTPtr->kmerArySI = 0;
   alnSTPtr->matrixSTPtr = 0;
   alnSTPtr->chainsSTPtr = 0;

   alnSTPtr->sizeSI = 0;
   blank_aln_mapRead(alnSTPtr);
} /*init_aln_mapRead*/

/*-------------------------------------------------------\
| Fun24: freeStack_aln_mapRead
|   - frees variables in an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer with variables to free
| Output:
|   - Modifies:
|     o frees kmerArySI, matrixSTPtr, and chainsSTPtr
|       in alnSTPtr
|     o sizeSI and lenSI in alnSTPtr to be 0
\-------------------------------------------------------*/
void
freeStack_aln_mapRead(
   struct aln_mapRead *alnSTPtr
){
   if(! alnSTPtr)
      return;

   if(alnSTPtr->kmerArySI)
      free(alnSTPtr->kmerArySI);
   if(alnSTPtr->matrixSTPtr)
      freeHeap_dirMatrix(alnSTPtr->matrixSTPtr);
   if(alnSTPtr->chainsSTPtr)
      freeHeap_chains_mapRead(alnSTPtr->chainsSTPtr);

   init_aln_mapRead(alnSTPtr);
} /*freeStack_aln_mapRead*/

/*-------------------------------------------------------\
| Fun25: freeHeap_aln_mapRead
|   - frees an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o frees alnSTPtr (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_aln_mapRead(
   struct aln_mapRead *alnSTPtr
){
   if(! alnSTPtr)
      return;

   freeStack_aln_mapRead(alnSTPtr);
   free(alnSTPtr);
} /*freeHeap_aln_mapRead*/

/*-------------------------------------------------------\
| Fun26: setup_aln_mapRead
|   - setups (memory allocate) an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer to setup
| Output:
|   - Modifies:
|     o allocates memory for kmerArySI, matrixSTPtr, and
|       chainsSTPtr in alnSTPtr
|     o sizeSI in alnSTPtr to be size of memory
|     o lenSI in alnSTPtr to be 0
|   - Returns:
|     o 0 for no errors (or no structure input)
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_aln_mapRead(
   struct aln_mapRead *alnSTPtr
){
   if(! alnSTPtr)
      return 0;

   freeStack_aln_mapRead(alnSTPtr);


   alnSTPtr->kmerArySI = malloc(256 * sizeof(signed int));
   if(! alnSTPtr->kmerArySI)
      goto memErr_fun23;
   alnSTPtr->sizeSI = 256;


   alnSTPtr->chainsSTPtr =
      mk_chains_mapRead(
         256, /*number chains, resize later if needed*/
         4096 /*will be resized later if needed*/
      ); /*calling this results in the blanking step
         `  being repeated twice, however,
         `   setup_aln_mapRead should only be called once
         */
   if(! alnSTPtr->chainsSTPtr)
      goto memErr_fun23;


   alnSTPtr->matrixSTPtr =
      malloc(sizeof(struct dirMatrix));
   if(! alnSTPtr->matrixSTPtr)
      goto memErr_fun23;
   init_dirMatrix(alnSTPtr->matrixSTPtr);


   blank_aln_mapRead(alnSTPtr);


   return 0;

   memErr_fun23:;
      return def_memErr_mapRead;
} /*setup_aln_mapRead*/

/*-------------------------------------------------------\
| Fun27: mk_aln_mapRead
|   - makes and returns heap allocated aln_mapRead struct
| Input:
| Output:
|   - Returns:
|     o pointer to new initialized and setup aln_mapRead
|       struct
|     o 0 for memory errors
\-------------------------------------------------------*/
struct aln_mapRead *
mk_aln_mapRead(
   void
){
   struct aln_mapRead *retSTPtr = 0;

   retSTPtr = malloc(sizeof(struct aln_mapRead));
   if(! retSTPtr)
      goto memErr_fun27;

   init_aln_mapRead(retSTPtr);
   if( setup_aln_mapRead(retSTPtr) )
      goto memErr_fun27;

   goto ret_fun27;

   memErr_fun27:;
      if(retSTPtr)
         freeHeap_aln_mapRead(retSTPtr);
      retSTPtr = 0;
      goto ret_fun27;

   ret_fun27:;
      return retSTPtr;
} /*mk_aln_mapRead*/

/*-------------------------------------------------------\
| Fun28: sortKmerIndex_mapRead
|   - sorts kmer index array by a kmer array
| Input:
|   - refKmerArySI:
|     o signed int array of reference kmers to use for
|       sorting
|   - indexArySI:
|     o signed int array to have sorted kmer index's
|     o needs to be same length as lenSI
|   - lenSI:
|     o length of kmer array (index 1)
| Output:
|   - Modifies:
|     o indexArySI to have index's sorted by refKmerArySI
|       - also sets index's for indexArySI
\-------------------------------------------------------*/
void
sortKmerIndex_mapRead(
   signed int *kmerArySI, /*kmers to make sorted index's*/
   signed int *indexArySI,/*get kmerArySI sorted index's*/
   signed int lenSI       /*array length (index 1)*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun28 TOC:
   '   - sorts kmer index array by a kmer array
   '   o fun28 sec01:
   '     - variable declerations
   '   o fun28 sec02:
   '     - find the number of rounds to sort for
   '   o fun28 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun28 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Number of sorting rounds*/
   signed int subSI = 0;
   signed int nextSI = 0;
   signed int lastSI = 0;
   signed int onSI = 0;

   /*Variables to incurment loops*/
   signed int siIndex = 0;
   signed int siElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun28 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   --lenSI; /*convert index 1 to index 0*/

   indexArySI[0] = 0;

   if(lenSI <= 0)
      return;

   /*set up index's*/
   for(
      subSI = 1;
      subSI < lenSI;
      ++subSI
   ) indexArySI[subSI] = (signed int) subSI;


   /*Recursion formsia: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subSI = 1; /*Initialzie first array*/

   while(subSI < lenSI)
      subSI = (3 * subSI) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun28 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subSI > 0)
   { /*Loop: all rounds*/
      for(
         siIndex = 0;
         siIndex <= subSI;
         ++siIndex
      ){ /*Loop: though sub array*/
         siElm = siIndex;

         for(
            siElm = siIndex;
            siElm + subSI <= lenSI;
            siElm += subSI
         ){ /*Loop: swap elements in subarray*/
            nextSI = siElm + subSI;

            if(
                 (unsigned int)
                    kmerArySI[ indexArySI[siElm] ]
               > (unsigned int)
                    kmerArySI[ indexArySI[nextSI] ]
            ) goto first_swap_fun28_sec03;

            else if(
                 (unsigned int)
                    kmerArySI[ indexArySI[siElm] ]
               < (unsigned int)
                    kmerArySI[ indexArySI[nextSI] ]
            ) ;

            else if(
               indexArySI[siElm] > indexArySI[nextSI]
            ){ /*Else If I need to swap an element*/
               first_swap_fun28_sec03:;

               indexArySI[siElm] ^= indexArySI[nextSI];
               indexArySI[nextSI] ^= indexArySI[siElm];
               indexArySI[siElm] ^= indexArySI[nextSI];

               lastSI = siElm;
               onSI = siElm;

               while(lastSI >= subSI)
               { /*loop: move swapped element back*/
                  lastSI -= subSI;

                  if(
                       (unsigned int)
                          kmerArySI[ indexArySI[onSI] ]
                     > (unsigned int)
                          kmerArySI[ indexArySI[lastSI] ]
                  ) break; /*Positioned the element*/

                  else if(
                       (unsigned int)
                          kmerArySI[ indexArySI[onSI] ]
                     < (unsigned int)
                          kmerArySI[ indexArySI[lastSI] ]
                  ) ; /*need to swap*/

                  else if(
                     indexArySI[onSI] > indexArySI[lastSI]
                  ) break; /*positioned element*/
   
                  indexArySI[onSI] ^= indexArySI[lastSI];
                  indexArySI[lastSI] ^= indexArySI[onSI];
                  indexArySI[onSI] ^= indexArySI[lastSI];

                  onSI = lastSI;
               } /*Loop: move swapped element back*/
            } /*Else If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subSI = (subSI - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*sortKmerIndex_mapRead*/

/*-------------------------------------------------------\
| Fun29: findKmer_mapRead
|   - finds first unsigned int in a range
| Input:
|   - indexArySI:
|     o signed int array of sorted index's to search
|     o make sure is sorted by sortKmerIndex_mapRead
|   - refKmerArySI:
|     o signed int array of reference kmers to search
|   - qrySI:
|     o query to find
|   - startSI:
|     o first kmer to start search at
|   - lenSI:
|     o length of array (index 1)
| Output:
|   - Returns:
|     o index of first item in array
|     o -1 if item not in array
\-------------------------------------------------------*/
signed long
findKmer_mapRead(
   signed int *indexArySI,   /*ref kmer index's*/
   signed int *refKmerArySI, /*reference kmers*/
   signed int qrySI,         /*query to search for*/
   signed long startSL,      /*first kmer to start at*/
   signed long lenSL         /*array length (index 1)*/
){
   signed long midSL = 0;

   --lenSL;

   if(lenSL <= 0)
      return -1;


   while(startSL <= lenSL)
   { /*Loop: find query*/
      midSL = (startSL + lenSL) >> 1;

      if(
          ( (unsigned int) qrySI )
        > ((unsigned int) refKmerArySI[indexArySI[midSL]])
      ) startSL = midSL + 1;

      else if(
          ( (unsigned int) qrySI )
        < ((unsigned int) refKmerArySI[indexArySI[midSL]])
      ) lenSL = midSL - 1;

      else
      { /*Else: is match*/
         if(! midSL)
            return midSL; /*at start of array*/

         if(
               ( (unsigned int) qrySI )
            == (  (unsigned int)
                  refKmerArySI[ indexArySI[midSL - 1] ]
               )
         ) lenSL = midSL - 1; /*still more matches*/

         else
            return midSL;
      } /*Else: is match*/
   } /*Loop: find query*/

   return -1; /*no matches found*/
} /*findKmer_mapRead*/

/*--------------------------------------------------------\
| Fun30: findChain_chains_mapRead
|   - finds chain with input kmer coordiantes in a
|     chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to chains_mapRead struct to search
|   - kmerUC:
|     o length of one kmer
|   - refSI:
|     o reference starting position of kmer
|   - qrySI:
|     o query starting position of kmer
| Output:
|   - Returns:
|     o index of at least one kmer chain with both
|       refSI and qrySI kmer coordaintes
|     o -1 if did not find kmer
\--------------------------------------------------------*/
signed long
findChain_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*has chains*/
   unsigned char kmerUC,    /*length of one kmer*/
   signed int refSI,        /*reference kmer start*/
   signed int qrySI         /*query kmer start*/
){
   signed long indexSL = 0;
   signed long refEndSL = 0;
   signed long qryEndSL = 0;

   if(chainsSTPtr->lenChainSI <= 0)
      return -1;

   /*originally I used a binary search, but then I
   `  realized that my data was not properly sorted
   `  due to an earlier kmer (ex position 130) consuming
   `  a later kmer (ex 203)
   */
   while(indexSL < chainsSTPtr->lenChainSI)
   { /*Loop: find kmers*/
      refEndSL = chainsSTPtr->refArySI[indexSL];
      refEndSL += chainsSTPtr->chainArySI[indexSL];
      refEndSL += kmerUC;
      --refEndSL;

      qryEndSL = chainsSTPtr->qryArySI[indexSL];
      qryEndSL += chainsSTPtr->chainArySI[indexSL];
      qryEndSL += kmerUC;
      --qryEndSL;

      if(qrySI >= qryEndSL)
         ++indexSL;
      else if(qrySI < chainsSTPtr->qryArySI[indexSL])
         ++indexSL;
      else if(refSI >= refEndSL)
         ++indexSL;
      else if(refSI < chainsSTPtr->refArySI[indexSL])
         ++indexSL;
      else
         return indexSL;
         /*chain includes reference and query kmers*/
   } /*Loop: find kmers*/

   return -1; /*nothing found*/
} /*findChain_chains_mapRead*/

/*-------------------------------------------------------\
| Fun31: kmerChain_mapRead
|   - merge kmers into chains (longer kmers)
| Input:
|   - chainsSTPtr:
|     o pointer to chainsSTPtr struct to add kmer chains
|       to
|   - minChainLenSI:
|     o minimum number of kmers to keep a chain
|   - kmerUC:
|     o length of one kmer
|   - qryKmerArySI:
|     o signed int array with query kmers (by position)
|   - lenQrySI:
|     o number kmers in qryKmerArySI
|   - refKmerArySI:
|     o signed int array with reference kmers to chain
|   - refIndexArySI:
|     o signed int array with refKmerArySI index's sorted
|       by kmer
|     o create with fun28 sortKmerIndex_mapRead
|   - lenRefSI:
|     o number kmers in refKmerArySI
|   - refStartSI:
|     o first refrence base in window
| Output:
|   - Modifies:
|     o refKmerArySI, qryKmerArySI, and chainLenArySI in
|       chainsSTPtr to have chain starting coordinates and
|       length
|     o scoreArySI and nextArySI in chainsSTPtr to have
|       all detected kmer positions filled with 0's
|     o indexArySI in chainsSTPtr to have index or -1 for
|       each query kmer start position in qryKmerArySI
|     o lenChainSI and lenIndexSI to have length of
|       indexArySI (num ref kmers) and number of chains
|     o sizeChainSI and sizeIndexSI in chainsSTPtr if
|       needed to rellocate memory for chain/index arrays
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|   - Note:
|     o sorting order of chains is query coordinate and
|       then reference coordinate
\-------------------------------------------------------*/
signed char
kmerChain_mapRead(
   struct chains_mapRead *chainsSTPtr,/*gets kmer chains*/
   signed int minChainLenSI,          /*min chain length*/
   unsigned int kmerUC,               /*length of 1 kmer*/
   signed int *qryKmerArySI,          /*query kmers*/
   signed int lenQrySI,               /*query length*/
   signed int *refKmerArySI,          /*reference kmers*/
   signed int *refIndexArySI,
     /*sorted refernce kmer index's (fun28)*/
   signed int lenRefSI,               /*number ref kmers*/
   signed int refStartSI              /*first ref base*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun31 TOC:
   '   - merge kmers to get longer kmers
   '   o fun31 sec01:
   '     - variable declarations
   '   o fun31 sec02:
   '     - initial memory allocation
   '   o fun31 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun31 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;          /*for returing errors*/

   signed long indexSL = 0;    /*index of kmer matches*/

   signed int qryIndexSI = 0;      /*qry kmer on*/
   signed int refIndexSI = 0;      /*table kmer on*/

   signed int lenChainSI = 0;      /*length of 1 chain*/
   signed int siKmer = 0;          /*reference kmer iter*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun31 Sec02:
   ^   - initial memory allocation
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_chains_mapRead(chainsSTPtr);

   if(chainsSTPtr->sizeIndexSI < lenQrySI)
   { /*If: need to resize index*/
      errSC =
         realloc_chains_mapRead(
            chainsSTPtr,
            0, /*not resizing chains array*/
            lenQrySI
         );
      if(errSC)
         goto memErr_fun31_sec04;
   } /*If: need to resize index*/

   chainsSTPtr->lenIndexSI = lenQrySI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun31 Sec03:
   ^   - find kmer chains
   ^   o fun31 sec03 sub01:
   ^     - find if have kmer chain(s) + start detect loop
   ^   o fun31 sec03 sub02:
   ^     - build all kmer chains for one ref kmer
   ^   o fun31 sec03 sub03:
   ^     - find next possible kmer match to ref kmer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun31 Sec03 Sub01:
   *   - find if have kmer chain(s) + start detect loop
   \*****************************************************/

   /* very hard to estimate O/theta (finding worst case
   `   is hard); this is an off the`cuff guess
   `  O[  q * log(r) qry-find * O[r * hits^2 * O[r|q
   `    + log(r) next-ref
   `   ]
   `   = O[q * r * log(2r + hits^2) * (r|q)]
   `      but r|q worst case is many duplicate kmers
   `      that can never meet the minimum chain
   `      length
   `   o no chain case, never extend chain beyond
   `     minimum chain length (mChain)
   `   = O[q * r * log(2r + hits^2) * mChain]
   `   o however, this assumes that for each q, we
   `     are iterating over the entireity of r, but
   `     this would require duplicate kmers, which
   `     would mean the minimum chain length,
   `     so there must be enough kmers that non-mathcing
   `     patterns form.
   `     * for 1 kmer only, it is dulpicate, everything
   `       is merged, so [r * log(hits^2) * mChain]
   `       becomes [1 * (1) + mChain]
   `     * for 2 kmers, query is divided by 2, inorder to
   `       be worst case the reference must have a
   `       different kmer at every mChain interval. This
   `       is highly unlikey in biological data. The max
   `       length before duplication is (2 ^ mChain)
   `     * howver, in in this case we have no hits, so
   `       log(2r + hits^2) becomes log(2r)
   `   = O[q * (r / 2) * log(2r + hits^2) * mChain]
   `   o the worst case is highly unlikely, also, it is
   `     very likely kmer diversity will be high, so this
   `     more likely scales closer, but not quite to
   `     kmer-space [r / min(length(r), 4 ^ kmer-length]
   `   = theta[q * kmer-space * log(2r + hits^2) * mChain]
   `   o for most cases the kmer space will be larger
   `     (7mer to 11mer; 16000 kmers to 4000000 kmers)
   `     than the read length, so only duplications
   `     will have an impact on the theta case, so
   `     log(2r + hits^2) might scale to log(r + hits^2)
   `   o the kmer scan step should reduce the chance
   `     of duplications, by limiting windows size
   `   o based on this it may often, very roughly scale to
   `   = theta[q * log(r + hits^2] * mChain]
   `     * probably my overally optimisitic guess
   `   o realisitically, hits will be few, so it the n^2
   `     step may be minimal, still it is nasty
   */

   for(
      siKmer = 0;
      siKmer < lenQrySI;
      ++siKmer
   ){ /*Loop: find matching kmers*/

      chainsSTPtr->indexArySI[siKmer] = -1;

      if(qryKmerArySI[siKmer] < 0)
      { /*If: anonymous kmer or end of kmers*/
         if(qryKmerArySI[siKmer] != def_endKmers_kmerFun)
            continue; /*anonymous base*/
         else
            break; /*no more query kmers*/
      } /*If: anonymous kmer or end of kmers*/

      indexSL =
         findKmer_mapRead(
            refIndexArySI,
            refKmerArySI,
            (unsigned int) qryKmerArySI[siKmer],
            0,
            lenRefSI
         ); /*find first query kmer for ref*/
            /*O[log(r)] ref-find*/
            /*unsigned int forces so negative values
            `  to be at end
            */

      if(indexSL < 0)
         continue; /*no match*/

      /**************************************************\
      * Fun31 Sec03 Sub02:
      *   - build all kmer chains for one ref kmer
      \**************************************************/

      /*  O[r] ref-loop
      ` * O[log(hits)] ref-find
      ` * O[r | q] extend
      */
      while(
            qryKmerArySI[siKmer]
         == refKmerArySI[ refIndexArySI[indexSL] ]
      ){ /*Loop: make kmer chains*/

         qryIndexSI = siKmer;
         refIndexSI = refIndexArySI[indexSL];
         lenChainSI = 0;

         if(
            findChain_chains_mapRead(
               chainsSTPtr,
               kmerUC,
               refIndexSI + refStartSI,
               qryIndexSI
            ) >= 0
         ) goto nextMatch_fun31_sec03_sub03;
           /*already have this kmer chain*/
           /*this is an n^2 step hit-find*/


         while(
                qryIndexSI < lenQrySI
             && refIndexSI < lenRefSI
             &&
                   qryKmerArySI[qryIndexSI++]
                == refKmerArySI[refIndexSI++]
         ) ++lenChainSI; /*find kmer chain length*/
           /*O[r | q] extend*/


         if(lenChainSI >= minChainLenSI)
         { /*If: keeping kmer chain*/

             if(
                   chainsSTPtr->lenChainSI
                >= chainsSTPtr->sizeChainSI
             ){ /*If: need more memory*/
                errSC =
                   realloc_chains_mapRead(
                      chainsSTPtr,
                      chainsSTPtr->lenChainSI << 1,
                      0 /*not rellocating index arrray*/
                   );
                if(errSC)
                   goto memErr_fun31_sec04;
             } /*If: need more memory*/

            if(chainsSTPtr->indexArySI[siKmer] < 0)
            { /*If: first chain for this query base*/
               chainsSTPtr->indexArySI[siKmer] =
                  chainsSTPtr->lenChainSI;

               chainsSTPtr->lastIndexSI = siKmer;
            } /*If: first chain for this query base*/

             /*reference kmer starting point*/
             chainsSTPtr->refArySI[
                chainsSTPtr->lenChainSI
             ] =  refIndexArySI[indexSL] + refStartSI;

             /*query starting point*/
             chainsSTPtr->qryArySI[
                chainsSTPtr->lenChainSI
             ] = siKmer;

             /*kmer chain length*/
             chainsSTPtr->chainArySI[
                chainsSTPtr->lenChainSI++
             ] = lenChainSI;
               
         } /*If: keeping kmer chain*/

         /***********************************************\
         * Fun31 Sec03 Sub03:
         *   - find next possible kmer match to ref kmer
         \***********************************************/

         nextMatch_fun31_sec03_sub03:;
            ++indexSL; /*get off last match*/

            if(indexSL >= lenRefSI)
               break; /*end of reference sequence*/
      } /*Loop: make kmer chains*/
   }  /*Loop: find matching kmers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun31 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun31_sec04;

   memErr_fun31_sec04:;
      errSC = def_memErr_mapRead;
      goto ret_fun31_sec04;

   ret_fun31_sec04:
      return errSC;
} /*kmerChain_mapRead*/

/*-------------------------------------------------------\
| Fun32: numQryChains_mapRead
|   - gets the number of chains assigned to a single
|     query kmer
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct with reference
|       chaines
|   - index:
|     o index of reference kmer (kmer number) to find
|       number of chains for
| Output:
|   - Returns:
|     o number of chains for reference kmer
\-------------------------------------------------------*/
#define numQryChains_mapRead(chainsSTPtr, index) ( ((chainsSTPtr)->indexArySI[(indexSL)] < 0) ? -1 : ( (index) < (chainsSTPtr)->lenIndexSI ) ? ( (chainsSTPtr)->indexArySI[(index) + 1] - (chainsSTPtr)->indexArySI[(index)] ) : ( (chainsSTPtr)->lenChainSI - (chainsSTPtr)->indexArySI[(index)] ) )
/*Logic:
`   - checks to see if reference kmer has any assigned
`     assigned chains
`   - turnary checks if have more index's (ref chains)
`     after the current index
`     (index < chainsSTPtr->lenIndexSI)
`   - if have more index's; length is found by
`     subtracing next reference starting index
`     (indexArySI[index + 1]) from the current index
`     (indexArySI[index])
`   - if at end; length is found by
`     subtracing the array length (chainsSTPtr->lenChainSI)
`     minus the references starting index
`     (chainsSTPtr->indexArySI[index])
*/

/*-------------------------------------------------------\
| Fun33: mergeChains_mapRead
|   - merges kmer chains from fun31 (kmerChain) into
|     the chain the maximizes base count
|   - not optimal solution, but also reduce complexity
|     by not counting for gaps
| Input:
|   - chainsSTPtr:
|     o chains_mapRead struct pionter with chains to
|       merge and to find "best" chain in
|   - lenKmerUC:
|     o size of one kmer
|   - gapSI:
|     o gap penalty to apply for indels between chains
|   - matchSI:
|     o score for match
|   - maxLenSI:
|     o maximum length between chains before merging
| Output:
|   - Modifies:
|     o nextArySI in chains_mapRead to have best scoring
|       trail for each chain
|     o scoreArySI to have most bases possible in a chain
|   - Returns:
|     o index of best merged chain
|     o def_noChains_mapRead if no chains in array
\-------------------------------------------------------*/ 
signed int
mergeChains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*chains to merge*/
   unsigned char lenKmerUC,            /*size of 1 kmer*/
   signed char gapSI,                  /*score for gaps*/
   signed char matchSI,                /*score for match*/
   signed int maxLenSI           /*length between chains*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun33 TOC:
   '   - merges kmer chains from fun31 (kmerChain)
   '   o fun33 sec01:
   '     - variable declarations
   '   o fun33 Sec02:
   '     - memory allocation
   '   o fun33 sec03:
   '     - check if one or no chains (mark and return)
   '   o fun33 sec04:
   '     - merge chains
   '   o fun33 sec05:
   '     - return results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun33 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int siChain = 0;    /*main chain on*/
   signed int siSub = 0;      /*subchain on*/
   signed int lastKmerSI = 0; /*index of last ref kmer*/
   signed int lastIndexSI = 0;/*index of last refkmer*/
   signed int subIndexSI = 0; /*index of subchain ref*/

   signed int scoreSI = 0;    /*score of merge*/
   signed int gapScoreSI = 0;
   signed int lenSI = 0;      /*number bases to add in*/

   signed int refUniqSI = 0;  /*number uniq bases in ref*/
   signed int refEndSI = 0;   /*ref ending coordinate*/

   signed int qryUniqSI = 0; /*number unique query bases*/
   signed int qryEndSI = 0;  /*query ending coordinate*/

   signed int maxNtSI = 0;    /*most bases in chain*/
   signed int curMaxNtSI = 0; /*current best base count*/
   signed int maxIndexSI = 0; /*index of chain maximized*/
   signed int curMaxIndexSI = 0;

   signed int secMaxNtSI = 0; /*2nd most bases in chain*/
   signed int curSecMaxNtSI = 0; /*current rounds score*/
      /*used as a quick check to reduces chances of a full
      `  scan. The logic is if
      `     (nt in chain) + secMaxNtSI
      `   < maxNtSI + (nt in chain) - overlap
      `  then I can safely say that I can not maximize
      `  the score any more for overlap cases
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun33 Sec02:
   ^   - memory allocation
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun33 Sec03:
   ^   - check if one or no chains (just mark and return)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(chainsSTPtr->lenChainSI <= 0)
      goto noChains_fun33_sec05;

   if(chainsSTPtr->lenChainSI == 1)
   { /*If: only one chain*/
      maxIndexSI = 0;

      chainsSTPtr->scoreArySI[0] =
         chainsSTPtr->chainArySI[0];
      chainsSTPtr->scoreArySI[0] += lenKmerUC;
      chainsSTPtr->scoreArySI[0] -= 1;

      chainsSTPtr->nextArySI[0] = 0; /*nothing after*/

      goto ret_fun33_sec05;
   } /*If: only one chain*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun33 Sec04:
   ^   - merge chains
   ^   o fun33 sec04 sub01:
   ^     - find last mapped reference kmer
   ^   o fun33 sec04 sub02:
   ^     - merge kmers
   ^   o fun33 sec04 sub03:
   ^     - score last reference kmer chains
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun33 Sec04 Sub01:
   *   - find last mapped reference kmer
   \*****************************************************/

   lastKmerSI = chainsSTPtr->lastIndexSI;

   if(lastKmerSI < 0)
      goto noChains_fun33_sec05;
   else
   { /*Else: have chains to merge*/
      lastIndexSI = lastKmerSI;
      lastKmerSI = chainsSTPtr->indexArySI[lastKmerSI];
   } /*Else: have chains to merge*/

   /*****************************************************\
   * Fun33 Sec04 Sub02:
   *   - score last reference kmer chains
   \*****************************************************/

   siChain = chainsSTPtr->lenChainSI - 1;

   while(
         chainsSTPtr->qryArySI[siChain]
      >= chainsSTPtr->qryArySI[lastKmerSI]
   ){ /*Loop: merge kmer chains*/
      chainsSTPtr->scoreArySI[siChain] =
         chainsSTPtr->chainArySI[siChain];

      chainsSTPtr->scoreArySI[siChain] += lenKmerUC;
      chainsSTPtr->scoreArySI[siChain] -= 1;
         /*account for index 1 of length (chainArySI)
         `   and position at first kmer nt
         */

      chainsSTPtr->nextArySI[siChain] = 0;/*end of merge*/

      /*want best score on first reference kmer, that
      `  way I avoid long gaps (hopefully)
      */
      if(chainsSTPtr->chainArySI[siChain] >= maxNtSI)
      { /*If: new best score*/
         secMaxNtSI = maxNtSI;
         maxNtSI = chainsSTPtr->scoreArySI[siChain];
         maxIndexSI = siChain;   
      } /*If: new best score*/

      else if(
         chainsSTPtr->scoreArySI[siChain] > secMaxNtSI
      ) secMaxNtSI = chainsSTPtr->scoreArySI[siChain];
         /*so always have second best case*/

      --siChain;

      if(siChain < 0)
         goto ret_fun33_sec05;
         /*no chains can be merged*/
   } /*Loop: merge kmer chains*/

   /*****************************************************\
   * Fun33 Sec04 Sub03:
   *   - score chains
   *   o fun33 sec04 sub03 cat01:
   *     - loop though chains & find end coords of chain
   *   o fun33 sec04 sub03 cat02:
   *     - handel no overlap cases (use highest score)
   *   o fun33 sec04 sub03 cat03:
   *     - handel overlap, but still best score cases
   *   o fun33 sec04 sub03 cat04:
   *     - handle cases were need to find best score
   *   o fun33 sec04 sub03 cat05:
   *     - merge valid subchains
   *   o fun33 sec04 sub03 cat06:
   *     - move to next sub chain
   *   o fun33 sec04 sub03 cat07:
   *     - make sure found best score is better than none
   *   o fun33 sec04 sub03 cat08:
   *     - move back one chain
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun33 Sec04 Sub03 Cat01:
   +   - loop though chains & find end coords of chain
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   curMaxNtSI = maxNtSI;
   curSecMaxNtSI = secMaxNtSI;
   curMaxIndexSI = maxIndexSI;

   if(! curSecMaxNtSI)
   { /*If: no secondary maximum chain*/
      curSecMaxNtSI = maxNtSI;
      secMaxNtSI = maxNtSI;
   } /*If: no secondary maximum chain*/

   /*this is kinda like the pheramon trail for the ant
   `   traveling salesman problem.
   */
   while(siChain >= 0)
   { /*Loop: merge chains*/
      lenSI = chainsSTPtr->chainArySI[siChain];
      lenSI += lenKmerUC;
      lenSI -= 1;

      qryEndSI = chainsSTPtr->qryArySI[siChain];
      qryEndSI += lenSI;
         /*account for index 1 of length (chainArySI)
         `   and position at first kmer nt
         */

      refEndSI = chainsSTPtr->refArySI[siChain];
      refEndSI += lenSI;
         /*account for index 1 of length (chainArySI)
         `   and position at first kmer nt
         */

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun33 Sec04 Sub03 Cat02:
      +   - handel no overlap cases (use highest score)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      qryUniqSI = chainsSTPtr->qryArySI[maxIndexSI];
      qryUniqSI -= qryEndSI;
      qryUniqSI *= gapSI;

      gapScoreSI = chainsSTPtr->refArySI[maxIndexSI];
      gapScoreSI -= refEndSI;
      gapScoreSI *= gapSI;

      gapScoreSI = min_genMath(qryUniqSI, gapScoreSI);

      qryUniqSI = qryEndSI;
      qryUniqSI -= chainsSTPtr->qryArySI[maxIndexSI];

      refUniqSI = refEndSI;
      refUniqSI -= chainsSTPtr->refArySI[maxIndexSI];

      lenSI *= matchSI;

      if(gapScoreSI + lenSI <= 0)
         goto findMerge_sec04_sub03_cat04;
         /*added bases do not cover added indels*/

      else if(
            qryUniqSI < 0
         && refUniqSI < 0
      ){ /*Else If: no overlap*/

         qryUniqSI *= -1;
         refUniqSI *= -1;

         /*make sure no large gaps*/
         if(qryUniqSI > maxLenSI )
            goto findMerge_sec04_sub03_cat04;
         else if(refUniqSI > maxLenSI )
            goto findMerge_sec04_sub03_cat04;

         chainsSTPtr->scoreArySI[siChain] =
            chainsSTPtr->chainArySI[siChain];
         chainsSTPtr->scoreArySI[siChain] += lenKmerUC;
         chainsSTPtr->scoreArySI[siChain] -= 1;
         chainsSTPtr->scoreArySI[siChain] += maxNtSI;

         chainsSTPtr->gapScoreArySI[siChain] =
            chainsSTPtr->gapScoreArySI[maxIndexSI];
         chainsSTPtr->gapScoreArySI[siChain] +=
            gapScoreSI;

         chainsSTPtr->nextArySI[siChain] = maxIndexSI;

         if(
               chainsSTPtr->scoreArySI[siChain]
            >= curMaxNtSI
         ){ /*If: new best score*/
            curSecMaxNtSI = curMaxNtSI;
            curMaxNtSI = chainsSTPtr->scoreArySI[siChain];
            curMaxIndexSI = siChain;
         } /*If: new best score*/

         else if(
              chainsSTPtr->scoreArySI[siChain]
            > curSecMaxNtSI
         ) curSecMaxNtSI=chainsSTPtr->scoreArySI[siChain];

         goto mvBack_fun33_sec04_sub03_cat08;
      }  /*Else If: no overlap*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun33 Sec04 Sub03 Cat03:
      +   - handel overlap, but still best score cases
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
            chainsSTPtr->qryArySI[siChain] + qryUniqSI
         >= qryEndSI
      ) goto findMerge_sec04_sub03_cat04;
        /*complete overlap, find next score*/

      else if(
            chainsSTPtr->refArySI[siChain] + refUniqSI
         >= refEndSI
      ) goto findMerge_sec04_sub03_cat04;
        /*complete overlap, find next score*/

      else if(
               chainsSTPtr->qryArySI[siChain]
            <= chainsSTPtr->qryArySI[maxIndexSI]
         && qryEndSI > chainsSTPtr->qryArySI[maxIndexSI]
      ) goto findMerge_sec04_sub03_cat04;
         /*complete overlap case*/
      else if(
               chainsSTPtr->refArySI[siChain]
            <= chainsSTPtr->refArySI[maxIndexSI]
         && refEndSI > chainsSTPtr->refArySI[maxIndexSI]
      ) goto findMerge_sec04_sub03_cat04;
         /*complete overlap case*/


      /*make sure no large gaps*/
      else if( (qryUniqSI * -1) > maxLenSI )
         goto findMerge_sec04_sub03_cat04;
      else if( (refUniqSI * -1) > maxLenSI )
         goto findMerge_sec04_sub03_cat04;

      else if(
              chainsSTPtr->chainArySI[siChain]
           - 1
           + lenKmerUC
           +  maxNtSI
           -  max_genMath(refUniqSI, qryUniqSI)
         > secMaxNtSI
      ){ /*Else If: with overlap still is best score*/
         chainsSTPtr->scoreArySI[siChain] =
            chainsSTPtr->chainArySI[siChain];
         chainsSTPtr->scoreArySI[siChain] -=
            max_genMath(
               refUniqSI,
               qryUniqSI
         ); /*subtract overlapping bases*/
         chainsSTPtr->scoreArySI[siChain] += maxNtSI;
         chainsSTPtr->scoreArySI[siChain] += lenKmerUC;
         chainsSTPtr->scoreArySI[siChain] -= 1;

         /*find number bases in chain*/
         qryEndSI -= chainsSTPtr->qryArySI[siChain];

         if(chainsSTPtr->scoreArySI[siChain] <= qryEndSI)
         { /*If: get better score without overlap*/
            chainsSTPtr->nextArySI[siChain] = 0;
            chainsSTPtr->scoreArySI[siChain] = qryEndSI;
         }  /*If: get better score without overlap*/

         else
            chainsSTPtr->nextArySI[siChain] = maxIndexSI;

         if(
              chainsSTPtr->scoreArySI[siChain]
            >= curMaxNtSI
         ){ /*If: new best score*/
            curSecMaxNtSI = curMaxNtSI;
            curMaxNtSI = chainsSTPtr->scoreArySI[siChain];
            curMaxIndexSI = siChain;
         } /*If: new best score*/

         else if(
              chainsSTPtr->scoreArySI[siChain]
            > curSecMaxNtSI
         ) curSecMaxNtSI=chainsSTPtr->scoreArySI[siChain];

         goto mvBack_fun33_sec04_sub03_cat08;
      }  /*Else If: with overlap still is best score*/
         
      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun33 Sec04 Sub03 Cat04:
      +   - handle cases were need to find best score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*this is the n^2 part, however, I am reallying
      `  on kmer space to be scarce enough that this step
      `  will be rare
      */

      findMerge_sec04_sub03_cat04:;

      siSub = lastKmerSI; /*last chain was on*/
      subIndexSI = lastIndexSI;/*qry index of last chain*/

      while(siSub < chainsSTPtr->lenChainSI)
      { /*Loop: scan next query kmer*/
         /*check if score improves alignment*/
         qryUniqSI = chainsSTPtr->qryArySI[siSub];
         qryUniqSI -= qryEndSI;
         qryUniqSI *= gapSI;

         gapScoreSI = chainsSTPtr->refArySI[siSub];
         gapScoreSI -= refEndSI;
         gapScoreSI *= gapSI;

         gapScoreSI =
            min_genMath(
               qryUniqSI,
               gapScoreSI
            );

         qryUniqSI = chainsSTPtr->qryArySI[siSub];
         qryUniqSI += chainsSTPtr->chainArySI[siSub];
         qryUniqSI += lenKmerUC;
         qryUniqSI -= 1;

         refUniqSI = chainsSTPtr->refArySI[siSub];
         refUniqSI += chainsSTPtr->chainArySI[siSub];
         refUniqSI += lenKmerUC;
         refUniqSI -= 1;


         /*skipping chains I can not merge*/
         if(refEndSI >= refUniqSI)
            ;
         else if(gapScoreSI + lenSI <= 0)
            goto nextSub_fun33_sec04_sub03_cat06;
         else if(qryEndSI >= qryUniqSI)
         { /*Else If: complete query overlap*/
             ++subIndexSI;

             if(subIndexSI >= chainsSTPtr->sizeIndexSI)
                break;

             while(
                   siSub
                >= chainsSTPtr->indexArySI[subIndexSI]
                /*siSub is always > -1, so works*/
             ){ /*Loop: find next open chain*/
                ++subIndexSI;

                if(subIndexSI >= chainsSTPtr->sizeIndexSI)
                   break;
             }  /*Loop: find next open chain*/

             if(subIndexSI >= chainsSTPtr->sizeIndexSI)
                break;

             siSub = chainsSTPtr->indexArySI[subIndexSI];
             continue;
         }  /*Else If: complete query overlap*/

         /*else if(
                  chainsSTPtr->qryArySI[siChain]
               <= chainsSTPtr->qryArySI[siSub]
            && qryEndSI > chainsSTPtr->qryArySI[siSub]
         ) ;

         else if(
                  chainsSTPtr->refArySI[siChain]
               <= chainsSTPtr->refArySI[siSub]
            && refEndSI > chainsSTPtr->refArySI[siSub]
         ) ;*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun33 Sec04 Sub03 Cat05:
         +   - merge valid subchains
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         else
         { /*Else: have gap or overlap (ref and query)*/
            qryUniqSI = qryEndSI;
            qryUniqSI -= chainsSTPtr->qryArySI[siSub];

            refUniqSI = refEndSI;
            refUniqSI -= chainsSTPtr->refArySI[siSub];

            scoreSI = chainsSTPtr->scoreArySI[siSub];
            scoreSI += chainsSTPtr->chainArySI[siChain];
            scoreSI += lenKmerUC;
            scoreSI -= 1;

            if(qryUniqSI * -1 > maxLenSI)
               goto nextSub_fun33_sec04_sub03_cat06;
            else if(refUniqSI * -1 > maxLenSI)
               goto nextSub_fun33_sec04_sub03_cat06;
            else if(
                  refUniqSI > 0
               || qryUniqSI > 0
            ){ /*If: have overlap*/
               scoreSI -= 
                  max_genMath( 
                     refUniqSI,
                     qryUniqSI
                  ); /*compinsate for overlap*/

               gapScoreSI -= 
                  max_genMath( 
                     refUniqSI,
                     qryUniqSI
                  ); /*compinsate for overlap*/

               if(gapScoreSI <= 0)
                  goto nextSub_fun33_sec04_sub03_cat06;
            }  /*If: have overlap*/

            else
            { /*Else: no overlap; remove large gaps*/
               qryUniqSI *= -1;
               refUniqSI *= -1;

               if(qryUniqSI > maxLenSI )
                  goto nextSub_fun33_sec04_sub03_cat06;
               else if(refUniqSI > maxLenSI )
                  goto nextSub_fun33_sec04_sub03_cat06;
            } /*Else: no overlap; remove large gaps*/

            if(scoreSI > chainsSTPtr->scoreArySI[siChain])
            { /*If: new local best score*/
               gapScoreSI += lenSI;
               chainsSTPtr->scoreArySI[siChain] = scoreSI;
               chainsSTPtr->nextArySI[siChain] = siSub;

               chainsSTPtr->gapScoreArySI[siChain] =
                  chainsSTPtr->gapScoreArySI[siSub];
               chainsSTPtr->gapScoreArySI[siChain] +=
                  gapScoreSI;

               if(scoreSI > curMaxNtSI)
               { /*If: new best score*/
                  curSecMaxNtSI = maxNtSI;
                  curMaxNtSI = scoreSI;
                  curMaxIndexSI = siChain;
               } /*If: new best score*/

               else if(scoreSI > curSecMaxNtSI)
                  curSecMaxNtSI = scoreSI;
            } /*If: new local best score*/

            else if(scoreSI > curSecMaxNtSI)
               curSecMaxNtSI = scoreSI;
         } /*Else: have gap or overlap (ref and query)*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun33 Sec04 Sub03 Cat06:
         +   - move to next sub chain
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         nextSub_fun33_sec04_sub03_cat06:;
            ++siSub; /*move to next hit*/

      } /*Loop: scan next query kmer*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun33 Sec04 Sub03 Cat07:
      +   - make sure found best score is better than none
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*find number bases in chain*/
      qryEndSI -= chainsSTPtr->qryArySI[siChain];

      if(! chainsSTPtr->scoreArySI[siChain])
      { /*If: no score not found*/
         chainsSTPtr->scoreArySI[siChain] =
            chainsSTPtr->chainArySI[siChain];
         chainsSTPtr->scoreArySI[siChain] += lenKmerUC;
         chainsSTPtr->scoreArySI[siChain] -= 1;

         chainsSTPtr->nextArySI[siChain] = 0;

         if(chainsSTPtr->scoreArySI[siChain] > curMaxNtSI)
         { /*If: new best score*/
            curSecMaxNtSI = maxNtSI;
            curMaxNtSI = chainsSTPtr->scoreArySI[siChain];
            curMaxIndexSI = siChain;
         } /*If: new best score*/

         else if(
              chainsSTPtr->scoreArySI[siChain]
            > curSecMaxNtSI
         ) curSecMaxNtSI = scoreSI;
      } /*If: no score not found*/

      else if(
         chainsSTPtr->scoreArySI[siChain] <= qryEndSI
      ){ /*If: get better score without overlap*/
         chainsSTPtr->nextArySI[siChain] = 0;
         chainsSTPtr->scoreArySI[siChain] = qryEndSI;
      }  /*If: get better score without overlap*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun33 Sec04 Sub03 Cat08:
      +   - move back one chain
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      mvBack_fun33_sec04_sub03_cat08:;
         --siChain;

         maxNtSI = curMaxNtSI;
         maxIndexSI = curMaxIndexSI;
         secMaxNtSI = curSecMaxNtSI;

         if(lastKmerSI <= 0)
            ; /*If: on last chain*/
         else if(siChain < 0)
            ; /*If: on last chain*/
         else if(
             chainsSTPtr->qryArySI[siChain]
           < chainsSTPtr->qryArySI[lastKmerSI - 1]
             /*lastKmerSI -1 lands on current kmer*/
         ){ /*If: on next query kmer*/
             /*update to add in new maximum values for the
             `  next reference kmer (round)
             */
            if(lastIndexSI > 0)
               --lastIndexSI;
            else
               goto ret_fun33_sec05; /*finished*/

            lastKmerSI =
               chainsSTPtr->indexArySI[lastIndexSI];

            while(lastKmerSI < 0)
            { /*Loop: find previous kmer*/
               if(lastIndexSI < 0)
                  goto ret_fun33_sec05; /*finished*/

               --lastIndexSI;

               lastKmerSI =
                  chainsSTPtr->indexArySI[lastIndexSI];
            } /*Loop: find previous kmer*/
         } /*If: on next query kmer*/
   } /*Loop: merge chains*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun33 Sec05:
   ^   - return results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   maxNtSI = curMaxNtSI;
   maxIndexSI = curMaxIndexSI;
   secMaxNtSI = curSecMaxNtSI;

   goto ret_fun33_sec05;

   noChains_fun33_sec05:;
      maxIndexSI = def_noChains_mapRead;
      goto ret_fun33_sec05;

   ret_fun33_sec05:;
      return maxIndexSI;
} /*mergeChains_mapRead*/

/*-------------------------------------------------------\
| Fun34: mergeToSam_mapRead
|   - convert output from mergeChains_map read to sam
|     file entry (mapq is set to 0)
| Input:
|   - chainsSTPtr:
|     o chains_mapRead struct pointer with chains to use
|       in making the alignment
|   - indexSI:
|     o starting index of merged chain to align
|   - alnEndsBl:
|     o align up to 1.5 chain lengths of the softmasked
|       ends of reads (0 = do not align)
|   - qrySTPtr:
|     o seqST struct with query sequence
|   - refSTPtr:
|     o seqST struct pointer with reference sequence
|   - lenKmerUC:
|     o length of one kmer in reference/query
|   - chainLenSI:
|     o minimum length of a chain
|     o used to align ends before chains, alignment size
|       is 1.5 chains (for chainLenSI = 20, is 30 bases)
|       * 1.5 is off the top of my head
|   - settings:
|     o settings for alignment
|   - matrixSTPtr:
|     o directional matrix to use in alignment step
|   - samSTPtr:
|     o samEntry struct pointer to add alignment to
| Output:
|   - Modifies:
|     o samSTPtr to have alignment; everything in samSTPtr
|       is modified
|     o matrixSTPtr (nothing usefull, here to reduce
|       memory allocations in alignment step)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
mergeToSam_mapRead(
   struct chains_mapRead *chainsSTPtr,/*chains to use*/
   signed int indexSI,         /*chain index to start*/
   signed char alnEndsBl,      /*1: map end masks*/
   struct seqST *qrySTPtr,     /*query sequence*/
   struct seqST *refSTPtr,     /*reference sequence*/
   unsigned char lenKmerUC,    /*length of one kmer*/
   signed int chainLenSI,      /*min length of one chain*/
   struct alnSet *settings,    /*aligment settings*/
   struct dirMatrix *matrixSTPtr,/*for alignment step*/
   struct samEntry *samSTPtr   /*gets final alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun34 TOC:
   '   - convert output form mergeChains_map read to sam
   '     file entry (mapq is set to 0)
   '   o fun34 sec01:
   '     - variable declarations
   '   o fun34 sec02:
   '     - initialize and setup
   '   o fun34 sec03:
   '     - align ending soft masked bases
   '   o fun34 sec04:
   '     - align ending soft masked bases
   '   o fun34 sec05:
   '     - finish up samEntry values (no mapq)
   '   o fun34 sec06:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun34 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   signed int lastChainBl = 0;
      /*tells if have done last chain yet*/

   /*coordiantes for aligning chunks*/
   signed int refStartSI = 0;
   signed int refEndSI = 0;
   signed int refLastEndSI = 0;

   signed int qryStartSI = 0;
   signed int qryEndSI = 0;
   signed int qryLastEndSI = 0;
   signed int qryLenSI = 0; /*for overlap indel detect*/
   signed int refLenSI = 0; /*for overlap indel detect*/

   /*for alignments*/
   unsigned int maskUI = 0; /*throw away*/
   unsigned int anonUI = 0; /*throw away*/
   unsigned int tmpUI = 0;  /*cigar lengths (throw away)*/

   /*for sam entry editing*/
   signed char *swapSCPtr = 0; /*for cigar reallocs*/
   signed int siCig = 0;
   signed int newCigPosSI = 0;

   /*for dealing with softmasking at ends*/
   signed int bestCigSI = 0;
   signed long scoreSL = 0;
   signed long bestScoreSL = 0;

   /*for dealing with soft masked bases at start*/
   signed int *cigHeapArySI = 0;
   signed char *cigHeapStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun34 Sec02:
   ^   - initialize and setup
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   chainLenSI += lenKmerUC; /*get bases in one chain*/
   chainLenSI += (chainLenSI >> 1); /*softmasks at end*/

   cigHeapArySI =
      malloc(
         ((chainLenSI << 1) + 8) * sizeof(signed int)
      );
   if(! cigHeapArySI)
      goto memErr_fun34_sec06;

   cigHeapStr =
      malloc(
         ((chainLenSI << 1) + 8) * sizeof(signed char)
      );
   if(! cigHeapStr)
      goto memErr_fun34_sec06;


   /*for merged chain gap filing*/
   blank_samEntry(samSTPtr);
   samSTPtr->flagUS = 0; /*mark forward map*/


   /*setup reference start coords*/
   samSTPtr->refStartUI =
      (unsigned int) chainsSTPtr->refArySI[indexSI];

   if(chainsSTPtr->qryArySI[indexSI] > 0)
   { /*If: have starting softmasking*/
      samSTPtr->cigTypeStr[0] = 'S';

      samSTPtr->cigArySI[0] =
         chainsSTPtr->qryArySI[indexSI];
         /*coordinates in index 0, but cigar index 1*/

      siCig = 1;
      samSTPtr->cigLenUI = 1;
   } /*If: have starting softmasking*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun34 Sec03:
   ^   - align gaps in merged chains
   ^   o fun34 sec03 sub01:
   ^     - setup for first chain in mapping
   ^   o fun34 sec03 sub02:
   ^     - expand cigar array memory + start loop
   ^   o fun34 sec03 sub03:
   ^     - find first and last base in current chain
   ^   o fun34 sec03 sub04:
   ^     - detect if have overlap in chains
   ^   o fun34 sec03 sub05:
   ^     - handel easy insertion/deletion cases
   ^   o fun34 sec03 sub06:
   ^     - handel easy (single) snp cases
   ^   o fun34 sec03 sub07:
   ^     - handel gaps that need to be aligned
   ^   o fun34 sec03 sub08:
   ^     - add chain length to matches
   ^   o fun34 sec03 sub09:
   ^     - move to next chain in merge
   ^   o fun34 sec03 sub10:
   ^     - check if have last chain to merge
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun34 Sec03 Sub01:
   *   - setup for first chain in mapping
   \*****************************************************/

   samSTPtr->cigTypeStr[siCig] = '=';
   samSTPtr->cigArySI[siCig] =
      chainsSTPtr->chainArySI[indexSI];
   samSTPtr->cigArySI[siCig] += (lenKmerUC - 1);
        /*last kmer is always two bases off*/
   ++samSTPtr->cigLenUI;
   ++siCig;

   qryLastEndSI = chainsSTPtr->qryArySI[indexSI];
   qryLastEndSI += chainsSTPtr->chainArySI[indexSI];
   qryLastEndSI += (lenKmerUC - 1);
     /*last kmer is always one base off (-1)
     `  and in index one (-1; total -2)
     */

   refLastEndSI = chainsSTPtr->refArySI[indexSI];
   refLastEndSI += chainsSTPtr->chainArySI[indexSI];
   refLastEndSI += (lenKmerUC - 1);
     /*last kmer is always one base off (-1)
     `  and in index one (-1; total -2)
     */

   if(chainsSTPtr->nextArySI[indexSI])
      goto nextChain_fun34_sec03_sub07;
         /*first index will often be 0*/
   else
   { /*Else: no next element*/
      samSTPtr->refEndUI = (unsigned int) refLastEndSI;
      samSTPtr->refStartUI =
         (unsigned int) chainsSTPtr->refArySI[indexSI];
      goto checkMask_fun34_sec03_sub09;
         /*only one good chain*/
   } /*Else: no next element*/

   /*****************************************************\
   * Fun34 Sec03 Sub02:
   *   - expand cigar array memory (if need) + start loop
   \*****************************************************/

   while(chainsSTPtr->nextArySI[indexSI] > 0)
   { /*Loop: while have chains to merge*/
      mergeChains_fun34_sec03_sub07:;

      /*doing + 3 because max for this loop is two
      `   entries, the thrid is for a softmask ending
      */
      if(
        siCig + 3 >= (signed int) samSTPtr->cigSizeUI)
      { /*If: need more cigar memory*/
         newCigPosSI= (signed int) samSTPtr->cigSizeUI;
         newCigPosSI += (newCigPosSI >> 1);

         swapSCPtr =
            realloc(
               samSTPtr->cigTypeStr,
               (newCigPosSI + 8) * sizeof(signed char)
            );
         if(! swapSCPtr)
            goto memErr_fun34_sec06;
         samSTPtr->cigTypeStr = swapSCPtr;

         swapSCPtr =
            (signed char *)
            realloc(
               samSTPtr->cigArySI,
               (newCigPosSI + 8) * sizeof(signed int)
            );
         if(! swapSCPtr)
            goto memErr_fun34_sec06;
         samSTPtr->cigArySI = (signed int *) swapSCPtr;

         samSTPtr->cigSizeUI =
            (unsigned int) newCigPosSI;
      } /*If: need more cigar memory*/

      /**************************************************\
      * Fun34 Sec03 Sub03:
      *   - find first and last base in current chain
      \**************************************************/

      /*find start and end of kmer chain*/
      qryStartSI = chainsSTPtr->qryArySI[indexSI];
      qryEndSI = qryStartSI;
      qryEndSI += chainsSTPtr->chainArySI[indexSI];
      qryEndSI += (lenKmerUC - 1);
        /*last kmer is always one base off (-1)
        `  and in index one (-1; total -2)
        */

      refStartSI = chainsSTPtr->refArySI[indexSI];
      refEndSI = refStartSI;
      refEndSI += chainsSTPtr->chainArySI[indexSI];
      refEndSI += (lenKmerUC - 1);
        /*last kmer is always one base off (-1)
        `  and in index one (-1; total -2)
        */

      /**************************************************\
      * Fun34 Sec03 Sub04:
      *   - detect if have overlap in chains
      \**************************************************/

      /*account for overlaps in kmer chains*/
      if(qryStartSI < qryLastEndSI)
      { /*If: have overlap*/
         chainOverlap_fun34_sec02_sub03:;

         qryLenSI = qryEndSI - qryLastEndSI;
         refLenSI = refEndSI - refLastEndSI;

         if(qryLenSI < refLenSI)
         { /*If: overlap caused a deletion*/
            samSTPtr->cigTypeStr[siCig] = 'D';
            samSTPtr->cigArySI[siCig] = refLenSI;
            samSTPtr->cigArySI[siCig] -= qryLenSI;

            ++siCig;
            qryStartSI = qryLastEndSI;
               /*account for overlap*/
         } /*If: overlap caused a deletion*/

         else if(qryLenSI > refLenSI)
         { /*If: overlap caused an indel*/
            samSTPtr->cigTypeStr[siCig] = 'I';
            samSTPtr->cigArySI[siCig] = qryLenSI;
            samSTPtr->cigArySI[siCig] -= refLenSI;

            ++siCig;
            refStartSI = refLastEndSI;
               /*account for overlap*/
         } /*If: overlap caused an deletion*/

         else
         { /*Else: overlap in bases (no indels)*/
            qryStartSI = qryLastEndSI + 1;
            refStartSI = qryLastEndSI + 1;
         } /*Else: overlap in bases (no indels)*/
         /*this case should never happen*/
      } /*If: have overlap*/

      else if(refStartSI < refLastEndSI)
         goto chainOverlap_fun34_sec02_sub03;
         /*means overlap between chains, this is most
         `  likely to happen for indels in homopolymers
         */

      /**************************************************\
      * Fun34 Sec03 Sub05:
      *   - handel easy insertion/deletion cases
      \**************************************************/

      /*check if have easy call*/
      else if((refStartSI - refLastEndSI) == 0)
      { /*If: no reference gap (means query gap)*/
         samSTPtr->cigTypeStr[siCig] = 'I';
         samSTPtr->cigArySI[siCig] =
            qryStartSI - qryLastEndSI; /*gap size*/

         ++siCig;
      } /*If: no reference gap (means query gap)*/

      else if((qryStartSI - qryLastEndSI) == 0)
      { /*If: no query gap (means reference gap)*/
         samSTPtr->cigTypeStr[siCig] = 'D';
         samSTPtr->cigArySI[siCig] =
            refStartSI - refLastEndSI; /*gap size*/

         ++siCig;
      } /*If: no query gap (means reference gap)*/

      /**************************************************\
      * Fun34 Sec03 Sub06:
      *   - handel easy (single) snp cases
      \**************************************************/

      else if(
            qryStartSI - qryLastEndSI == 1
         && refStartSI - refLastEndSI == 1
      ){ /*Else If: mismatch*/
         samSTPtr->cigTypeStr[siCig] = 'X';
         samSTPtr->cigArySI[siCig] = 1;
         ++siCig;
      }  /*Else If: mismatch*/

      /**************************************************\
      * Fun34 Sec03 Sub07:
      *   - handel gaps that need to be aligned
      \**************************************************/

      else
      { /*Else: need to align*/
         qrySTPtr->offsetSL = qryLastEndSI;
         qrySTPtr->endAlnSL = qryStartSI;
         --qrySTPtr->endAlnSL; /*avoid 1st base in chain*/

         refSTPtr->offsetSL = refLastEndSI;

         refSTPtr->endAlnSL = (unsigned long) refStartSI;
         --refSTPtr->endAlnSL; /*avoid 1st base in chain*/

         needle(
            qrySTPtr,
            refSTPtr,
            matrixSTPtr,
            settings
         );

         if(matrixSTPtr->errSC == def_memErr_needle)
            goto memErr_fun34_sec06;

         samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';
            /*so getCig_dirMatrix knows is new cigar*/

         newCigPosSI =
            getCig_dirMatrix(
              matrixSTPtr,           /*has alignment*/
              0,                     /*use matrix index*/
              0,                     /*always forward*/
              qrySTPtr,              /*query sequence*/
              refSTPtr,              /*reference sequence*/
              &samSTPtr->cigTypeStr, /*cigar type array*/
              &samSTPtr->cigArySI,  /*has cigar lengths*/
              samSTPtr->cigLenUI,   /*position in cigar*/
              &samSTPtr->cigSizeUI, /*length of cigar*/
              (unsigned int *) &refLastEndSI,/*aln start*/
              &anonUI,               /*ignoring*/
              &maskUI,               /*ignoring*/
              1,                     /*keep ending dels*/
              settings               /*settings for aln*/
            ); /*add aligment to cigar*/

         siCig = newCigPosSI + 1;

         if( newCigPosSI < 0)
            goto memErr_fun34_sec06;

      } /*Else: need to align*/

      /**************************************************\
      * Fun34 Sec03 Sub08:
      *   - add chain length to matches
      \**************************************************/

      samSTPtr->refEndUI = (unsigned int) refEndSI;
         /*mark end of alignment*/
      samSTPtr->readLenUI = (unsigned int) qryEndSI;
         /*so know last query base*/

      /*setup for next round*/
      qryLastEndSI = qryEndSI;
      refLastEndSI = refEndSI;

      /*add chain lengths in*/
      qryEndSI -= qryStartSI; /*find query length*/
      refEndSI -= refStartSI; /*find reference length*/
         /*already accounted for extra bases in last
         `  kmer when found qryEndSI and refEndSI
         */

      samSTPtr->cigArySI[siCig] =
         min_genMath(
            qryEndSI,
            refEndSI
         );

      samSTPtr->cigTypeStr[siCig++] = '=';
      samSTPtr->cigLenUI = (unsigned int) siCig;

      /**************************************************\
      * Fun34 Sec03 Sub09:
      *   - move to next chain in merge
      \**************************************************/

      nextChain_fun34_sec03_sub07:;
         if(chainsSTPtr->nextArySI[indexSI])
            indexSI = chainsSTPtr->nextArySI[indexSI];
            /*move to next chain in merge*/
   } /*Loop: while have chains to merge*/

   /*****************************************************\
   * Fun34 Sec03 Sub10:
   *   - check if have last chain to merge
   \*****************************************************/

   if(! lastChainBl)
   { /*If: still need to finish the last chain*/
      lastChainBl = 1;
      goto mergeChains_fun34_sec03_sub07;
   } /*If: still need to finish the last chain*/

   checkMask_fun34_sec03_sub09:;

   samSTPtr->cigTypeStr[siCig] = '\0';

   if(qryLastEndSI < qrySTPtr->seqLenSL)
      /*< ( (signed int) (qrySTPtr->seqLenSL - 1) )*/
   { /*If: have softmasking at end*/
      /*cig softmasking length at end*/
      samSTPtr->cigArySI[siCig] =
        (signed int) qrySTPtr->seqLenSL;

      samSTPtr->cigArySI[siCig] -= qryLastEndSI;

      samSTPtr->cigTypeStr[siCig++] = 'S';
      samSTPtr->cigTypeStr[siCig] = '\0';
      ++samSTPtr->cigLenUI;

      samSTPtr->readLenUI =
         (signed int) qrySTPtr->seqLenSL;
   } /*If: have softmasking at end*/

   if(! alnEndsBl)
      goto skipEndMask_fun34_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun34 Sec04:
   ^   - align starting soft masked bases
   ^   o fun34 sec04 sub01:
   ^     - deal with starting softmasking
   ^   o fun34 sec04 sub02:
   ^     - deal with ending softmasking
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun34 Sec04 Sub01:
   *   - deal with starting softmasking
   *   o fun34 sec04 sub01 cat01:
   *     - find softamsking coordinates to align
   *   o fun34 sec04 sub01 cat02:
   *     - do needleman alignment8
   *   o fun34 sec04 sub01 cat03:
   *     - get cigar and check if keeping
   *   o fun34 sec04 sub01 cat04:
   *     - check if keeping cigar
   *   o fun34 sec04 sub01 cat05:
   *     - realloc main memory (if more needed)
   *   o fun34 sec04 sub01 cat06:
   *     - make room in old cigar entry for extra entry
   *   o fun34 sec04 sub01 cat07:
   *     - add new cigar entries in
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun34 Sec04 Sub01 Cat01:
   +   - find softamsking coordinates to align
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(samSTPtr->cigTypeStr[0] == 'S')
   { /*If: have softmasking at start*/
      /*adding to end so can do swap*/

      if(samSTPtr->refStartUI < (unsigned int) chainLenSI)
         refStartSI = 0;
      else
         refStartSI =
            (signed int)
            samSTPtr->refStartUI - chainLenSI;

      if(samSTPtr->cigArySI[0] < chainLenSI)
         qryStartSI = 0;
      else
         qryStartSI =
            (signed int)
            samSTPtr->cigArySI[0] - chainLenSI;

      if(samSTPtr->refStartUI > 0)
         refEndSI = samSTPtr->refStartUI - 1;
      else
         refEndSI = 0;

      qryEndSI = samSTPtr->cigArySI[0];
      --qryEndSI;


      /*check if both sequences have softmasking*/
      if(refEndSI <= refStartSI)
         goto skipStartMask_fun34_sec04_sub02;
      else if (qryEndSI <= qryStartSI)
         goto skipStartMask_fun34_sec04_sub02;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub01 Cat02:
      +   - do needleman alignment
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*this is here to handel cases were reference is
      `   much longer than query
      */
      qryLastEndSI = qryEndSI - qryStartSI;
      refLastEndSI = refEndSI - refStartSI;

      bestScoreSL = qryLastEndSI;
      bestScoreSL += (qryLastEndSI >> 1);
         /*qry length + 50% query length (for indels)*/

      if(bestScoreSL < refLastEndSI)
         refStartSI = refEndSI - bestScoreSL;
         /*qry is much shorter than reference*/
      else
      { /*Else: need to check if query is longer*/
         bestScoreSL = refLastEndSI;
         bestScoreSL += (refLastEndSI >> 2);
            /*ref length + 25% query length (for indels)*/

         if(bestScoreSL < qryLastEndSI)
            qryStartSI = qryEndSI - bestScoreSL;
            /*ref is much shorter than query*/
      } /*Else: need to check if query is longer*/


      qrySTPtr->offsetSL = qryStartSI;
      qrySTPtr->endAlnSL = qryEndSI;

      refSTPtr->offsetSL = refStartSI;
      refSTPtr->endAlnSL = refEndSI;

       needle(
         qrySTPtr,
         refSTPtr,
         matrixSTPtr,
         settings
      );

      if(matrixSTPtr->errSC == def_memErr_needle)
         goto memErr_fun34_sec06;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub01 Cat03:
      +   - get cigar
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      refLastEndSI = refStartSI;
      tmpUI = (chainLenSI << 1);
      cigHeapStr[0] = '\0';

      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';
         /*so getCig_dirMatrix knows is new cigar*/

      newCigPosSI =
         getCig_dirMatrix(
           matrixSTPtr,           /*has alignment*/
           0,                     /*use matrix index*/
           0,                     /*always forward*/
           qrySTPtr,              /*query sequence*/
           refSTPtr,              /*reference sequence*/
           &cigHeapStr,           /*cigar type array*/
           &cigHeapArySI,         /*has cigar lengths*/
           0,                     /*first cigar pos (0)*/
           &tmpUI,                /*cigar buffer size*/
           (unsigned int *) &refLastEndSI,/*start of aln*/
           &anonUI,               /*ignoring*/
           &maskUI,               /*ignoring*/
           0,                     /*keep ending dels*/
           settings               /*settings for aln*/
         ); /*add aligment to cigar*/

      if(newCigPosSI < 0)
         goto memErr_fun34_sec06;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub01 Cat04:
      +   - check if keeping cigar
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      siCig = (unsigned int) newCigPosSI;
      scoreSL = 0;
      bestScoreSL = 0;
      bestCigSI = 0;

      refStartSI = refSTPtr->endAlnSL;
      qryStartSI = qrySTPtr->endAlnSL;

      /*needleman only reutrns match/snp/ins/del*/
      while(siCig >= 0)
      { /*Loop: find rough score*/
         if(cigHeapStr[siCig] == '=')
            goto matchSnpStart_fun34_sec04_sub01_cat03;

         else if(cigHeapStr[siCig] == 'X')
         { /*Else If: snp/match case*/
            matchSnpStart_fun34_sec04_sub01_cat03:;

            tmpUI = (unsigned int) cigHeapArySI[siCig];

            while(tmpUI > 0)
            { /*Loop: get score*/
               scoreSL +=
                  getScore_alnSet(
                     qrySTPtr->seqStr[qryStartSI],
                     refSTPtr->seqStr[refStartSI],
                     settings
                  );

               --qryStartSI;
               --refStartSI;
               --tmpUI;
            } /*Loop: get score*/

            if(scoreSL > bestScoreSL)
            { /*If: new best score*/
               bestScoreSL = scoreSL;
               bestCigSI = siCig;
            } /*If: new best score*/
         } /*Else If: snp/match case*/

         else if(cigHeapStr[siCig] == 'I')
         { /*Else If: insertion*/
            qryStartSI -= cigHeapArySI[siCig];

            goto indelStart_fun34_sec04_sub01_cat03;
         } /*Else If: insertion*/

         else
         { /*Else: deletion*/
            refStartSI -= cigHeapArySI[siCig];

            indelStart_fun34_sec04_sub01_cat03:;

            scoreSL += settings->gapSS;
            scoreSL +=
              (
                   settings->extendSS
                 * (cigHeapArySI[siCig] - 1)
              );
         } /*Else: deletion*/

         --siCig;
      } /*Loop: find rough score*/

      if(bestScoreSL <= 0)
         goto skipStartMask_fun34_sec04_sub02;

      ++newCigPosSI; /*convert to index 1*/
      maskUI = (unsigned int) (newCigPosSI - bestCigSI);

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub01 Cat05:
      +   - realloc main memory (if more need)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      indexSI = (signed int) samSTPtr->cigLenUI + maskUI;

      if(
            indexSI + 1
         >= (signed int) samSTPtr->cigSizeUI
      ){ /*If: need more memory*/

         swapSCPtr =
            realloc(
               samSTPtr->cigTypeStr,
               (indexSI + 8) * sizeof(signed char)
            );
         if(! swapSCPtr)
            goto memErr_fun34_sec06;
         samSTPtr->cigTypeStr = swapSCPtr;

         swapSCPtr =
            (signed char *)
            realloc(
               samSTPtr->cigArySI,
               (indexSI + 8) * sizeof(signed int)
            );
         if(! swapSCPtr)
            goto memErr_fun34_sec06;
         samSTPtr->cigArySI = (signed int *) swapSCPtr;

         samSTPtr->cigSizeUI = (unsigned int) indexSI;
      } /*If: need more memory*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub01 Cat06:
      +   - make room in old cigar entry for extra entry
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      indexSI = (signed int) samSTPtr->cigLenUI - 1;
      indexSI -= (! bestCigSI);
         /*reduce copy range by one if no softmasking*/
      indexSI += maskUI;
         /*removes softmasking at start if needed*/

      /*deal with softmasking not aligined*/
      if(bestCigSI == 1)
      { /*If: first entry could be a deletion*/
         if(qrySTPtr->offsetSL)
            bestCigSI = 1;

         else if(cigHeapStr[0] == 'D')
         { /*If: no query bases at first entry*/
            --indexSI;
            bestCigSI = 0;
         } /*If: no query bases at first entry*/
      } /*If: first entry could be a deletion*/

      else if(bestCigSI)
         ;
      else if(qrySTPtr->offsetSL)
      { /*Else If: did not align all query bases*/
         ++indexSI;
         bestCigSI = 1;
      } /*Else If: did not align all query bases*/

      siCig = samSTPtr->cigLenUI - 1;
         /*-1 to conver index 1 to index 0*/
      samSTPtr->cigLenUI = indexSI + 1;
      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';

      while(siCig > 0)
      { /*Loop: move cigar up*/
         samSTPtr->cigTypeStr[indexSI] =
            samSTPtr->cigTypeStr[siCig];
         samSTPtr->cigArySI[indexSI] =
            samSTPtr->cigArySI[siCig];

         --siCig;
         --indexSI;
      } /*Loop: move cigar up*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub01 Cat07:
      +   - add new cigar entries in
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      indexSI = newCigPosSI - (signed int) maskUI;
         /*find first cigar entry to copy*/
 
      siCig = 0;

      siCig = !! bestCigSI;
         /*goes to 1 if index of best sub-alignment
         `   (bestCigSI) is not 0 (sofmasked bases),
         `   else is 0 (got complete alignment)
         */

      maskUI = 0;

      while(indexSI < newCigPosSI)
      { /*Loop: add new cigar entries in*/
         samSTPtr->cigTypeStr[siCig]= cigHeapStr[indexSI];
         samSTPtr->cigArySI[siCig]=cigHeapArySI[indexSI];

         if(cigHeapStr[indexSI] != 'I')
            samSTPtr->refStartUI -= cigHeapArySI[indexSI];
            /*adjust reference starting position*/

         if(
               bestCigSI
            && cigHeapStr[indexSI] != 'D'
         ) samSTPtr->cigArySI[0] -= cigHeapArySI[indexSI];
           /*if need to adjust sofmask count*/

         ++siCig;
         ++indexSI;
      } /*Loop: add new cigar entries in*/
   } /*If: have softmasking at start*/

   /*****************************************************\
   * Fun34 Sec04 Sub02:
   *   - deal with ending softmasking
   *   o fun34 sec04 sub02 cat01:
   *     - find softamsking coordinates to align
   *   o fun34 sec04 sub02 cat02:
   *     - do needleman alignment8
   *   o fun34 sec04 sub02 cat03:
   *     - get cigar and check if keeping
   *   o fun34 sec04 sub02 cat04:
   *     - realloc main memory (if more needed)
   *   o fun34 sec04 sub02 cat05:
   *     - add new cigar to samEntry
   *   o fun34 sec04 sub02 cat06:
   *     - adjust softmasking at start
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun34 Sec04 Sub02 Cat01:
   +   - find softamsking coordinates to align
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   skipStartMask_fun34_sec04_sub02:;

   siCig = samSTPtr->cigLenUI - 1;

   if(samSTPtr->cigTypeStr[siCig] == 'S')
   { /*If: have ending softmasking*/
      /*adding to end so can do swap*/

      refStartSI = (signed int) samSTPtr->refEndUI;
      refEndSI = refStartSI;
      refEndSI += chainLenSI;
      --refEndSI; /*convert index one to index zero*/

      if(refEndSI >= (signed int) refSTPtr->seqLenSL)
         refEndSI = (signed int) refSTPtr->seqLenSL - 1;

      qryStartSI = (signed int) qrySTPtr->seqLenSL;
      qryStartSI -= samSTPtr->cigArySI[siCig];
      qryEndSI = qryStartSI;
      qryEndSI += chainLenSI;
      --qryEndSI; /*convert index 1 to index 0*/

      if(qryEndSI >= (signed int) qrySTPtr->seqLenSL)
         qryEndSI = (signed int) qrySTPtr->seqLenSL - 1;

      /*check if have sequence to align*/
      if(refEndSI <= refStartSI)
        goto skipEndMask_fun34_sec05;
      if(qryEndSI <= refStartSI)
        goto skipEndMask_fun34_sec05;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub02 Cat02:
      +   - do needleman alignment
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*this is here to handel cases were reference is
      `   much longer than query
      */
      qryLastEndSI = qryEndSI - qryStartSI;
      refLastEndSI = refEndSI - refStartSI;

      bestScoreSL = qryLastEndSI;
      bestScoreSL += (qryLastEndSI >> 1);
         /*qry length + 50% query length (for indels)*/

      if(bestScoreSL < refLastEndSI)
         refEndSI = refStartSI + bestScoreSL;
         /*qry is much shorter than reference*/
      else
      { /*Else: need to check if query is longer*/
         bestScoreSL = refLastEndSI;
         bestScoreSL += (refLastEndSI >> 2);
            /*ref length + 25% query length (for indels)*/

         if(bestScoreSL < qryLastEndSI)
            qryEndSI = qryStartSI + bestScoreSL;
            /*ref is much shorter than query*/
      } /*Else: need to check if query is longer*/


      qrySTPtr->offsetSL = qryStartSI;
      qrySTPtr->endAlnSL = qryEndSI;

      refSTPtr->offsetSL = refStartSI;
      refSTPtr->endAlnSL = refEndSI;

      needle(
         qrySTPtr,
         refSTPtr,
         matrixSTPtr,
         settings
      );

      if(matrixSTPtr->errSC == def_memErr_needle)
         goto memErr_fun34_sec06;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub02 Cat03:
      +   - get cigar and check if keeping
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      refLastEndSI = refStartSI;
      tmpUI = chainLenSI << 1;
      cigHeapStr[0] = '\0';

      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';
         /*so getCig_dirMatrix knows is new cigar*/

      newCigPosSI =
         getCig_dirMatrix(
            matrixSTPtr,         /*has alignment*/
            0,                   /*use matrix index*/
            0,                   /*always forward*/
            qrySTPtr,            /*query sequence*/
            refSTPtr,            /*reference sequence*/
            &cigHeapStr,         /*cigar type array*/
            &cigHeapArySI,       /*has cigar lengths*/
            0,                     /*first cigar pos (0)*/
            &tmpUI,                /*cigar buffer size*/
            (unsigned int *) &refStartSI,/*start of aln*/
            &anonUI,             /*ignoring*/
            &maskUI,             /*ignoring*/
            1,                   /*keep ending dels*/
            settings             /*settings for aln*/
         ); /*add aligment to cigar*/
      
      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub02 Cat04:
      +   - check if keeping cigar
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      siCig = 0;
      scoreSL = 0;
      bestScoreSL = 0;
      bestCigSI = 0;

      refStartSI = refSTPtr->offsetSL;
      qryStartSI = qrySTPtr->offsetSL;

      /*needleman only reutrns match/snp/ins/del*/
      while(siCig <= newCigPosSI)
      { /*Loop: find rough score*/
         if(cigHeapStr[siCig] == '=')
            goto matchSnpEnd_fun34_sec04_sub02_cat03;

         else if(cigHeapStr[siCig] == 'X')
         { /*Else If: snp/match case*/
            matchSnpEnd_fun34_sec04_sub02_cat03:;

            tmpUI = (unsigned int) cigHeapArySI[siCig];

            while(tmpUI > 0)
            { /*Loop: get score*/
               scoreSL +=
                  getScore_alnSet(
                     qrySTPtr->seqStr[qryStartSI],
                     refSTPtr->seqStr[refStartSI],
                     settings
                  );

               ++qryStartSI;
               ++refStartSI;
               --tmpUI;
            } /*Loop: get score*/

            if(scoreSL > bestScoreSL)
            { /*If: new best score*/
               bestScoreSL = scoreSL;
               bestCigSI = siCig;
            } /*If: new best score*/
         } /*Else If: snp/match case*/

         else if(cigHeapStr[siCig] == 'I')
         { /*Else If: insertion*/
            qryStartSI += cigHeapArySI[siCig];

            goto indelEnd_fun34_sec04_sub02_cat03;
         } /*Else If: insertion*/

         else
         { /*Else: deletion*/
            refStartSI += cigHeapArySI[siCig];

            indelEnd_fun34_sec04_sub02_cat03:;

            scoreSL += settings->gapSS;
            scoreSL +=
              (
                   settings->extendSS
                 * (cigHeapArySI[siCig] - 1)
              );
         } /*Else: deletion*/

         ++siCig;
      } /*Loop: find rough score*/

      if(bestScoreSL <= 0)
         goto skipEndMask_fun34_sec05;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub02 Cat05:
      +   - realloc main memory (if more needed)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      indexSI =
         (signed int) samSTPtr->cigLenUI + bestCigSI;

      if(
            indexSI + 1
         >= (signed int) samSTPtr->cigSizeUI
      ){ /*If: need more memory*/

         swapSCPtr =
            realloc(
               samSTPtr->cigTypeStr,
               (indexSI + 8) * sizeof(signed char)
            );
         if(! swapSCPtr)
            goto memErr_fun34_sec06;
         samSTPtr->cigTypeStr = swapSCPtr;

         swapSCPtr =
            (signed char *)
            realloc(
               samSTPtr->cigArySI,
               (indexSI + 8) * sizeof(signed int)
            );
         if(! swapSCPtr)
            goto memErr_fun34_sec06;

         samSTPtr->cigArySI = (signed int *) swapSCPtr;
         samSTPtr->cigSizeUI = (unsigned int) indexSI;
      } /*If: need more memory*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun34 Sec04 Sub02 Cat06:
      +   - add new cigar to samEntry
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      siCig = 0;
      --samSTPtr->cigLenUI;
         /*so will override softmasked base entry at end*/
      qryStartSI = qrySTPtr->offsetSL;

      while(siCig <= (signed int) bestCigSI)
      { /*Loop: move cigar up*/
         
         if(cigHeapStr[siCig] != 'D')
            qryStartSI += cigHeapArySI[siCig];
               /*allows me to find how many softmasked
               `  bases at end
               */

         if(cigHeapStr[siCig] != 'I')
            samSTPtr->refEndUI += cigHeapArySI[siCig];
               /*reference length*/

         samSTPtr->cigTypeStr[samSTPtr->cigLenUI] =
            cigHeapStr[siCig];

         samSTPtr->cigArySI[samSTPtr->cigLenUI] =
            cigHeapArySI[siCig];


         ++siCig;
         ++samSTPtr->cigLenUI;
      } /*Loop: move cigar up*/

      if(qryStartSI < (signed int) qrySTPtr->seqLenSL)
      { /*If: had softmasking*/
         samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = 'S';
         samSTPtr->cigArySI[samSTPtr->cigLenUI] = 
              (signed int) qrySTPtr->seqLenSL
            - qryStartSI;
         ++samSTPtr->cigLenUI;
      } /*If: had softmasking*/

      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';
   } /*If: have ending softmasking*/
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun34 Sec05:
   ^   - finish up samEntry values (no mapq)
   ^   o fun34 sec05 sub01:
   ^     - add read lengths
   ^   o fun34 sec05 sub02:
   ^     - copy read ids
   ^   o fun34 sec05 sub03:
   ^     - copy sequence
   ^   o fun34 sec05 sub04:
   ^     - copy q-score
   ^   o fun34 sec05 sub05:
   ^     - find mutation counts
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   skipEndMask_fun34_sec05:;

   /*****************************************************\
   * Fun34 Sec05 Sub01:
   *   - add read lengths
   \*****************************************************/

   /*add in read lengths*/
   samSTPtr->readLenUI =
      (unsigned int) qrySTPtr->seqLenSL;

   samSTPtr->alnReadLenUI =
      samSTPtr->refEndUI - samSTPtr->refStartUI;

   /*****************************************************\
   * Fun34 Sec05 Sub02:
   *   - copy read ids
   \*****************************************************/

   if(refSTPtr->idStr[0] == '@')
      indexSI = 1;
   else if(refSTPtr->idStr[0] == '>')
      indexSI = 1;
   else
      indexSI = 0;

   /*copy read ids*/
   samSTPtr->refIdLenUC =
      (unsigned int)
      cpWhite_ulCp(
         samSTPtr->refIdStr,
         &refSTPtr->idStr[indexSI]
      );


   if(qrySTPtr->idStr[0] == '@')
      indexSI = 1;
   else if(qrySTPtr->idStr[0] == '>')
      indexSI = 1;
   else
      indexSI = 0;

   samSTPtr->qryIdLenUC =
      (unsigned int)
      cpWhite_ulCp(
         samSTPtr->qryIdStr,
         &qrySTPtr->idStr[indexSI]
      );

   /*****************************************************\
   * Fun34 Sec05 Sub03:
   *   - copy sequence
   \*****************************************************/

   if(samSTPtr->seqSizeUI < qrySTPtr->seqLenSL)
   { /*If: need more memory for sequence*/
      free(samSTPtr->seqStr);
      samSTPtr->seqStr = 0;

      samSTPtr->seqStr =
         malloc(
              (qrySTPtr->seqLenSL + 9)
            * sizeof(unsigned char)
         );

      if(! samSTPtr)
         goto memErr_fun34_sec06;

      samSTPtr->seqSizeUI =
        (unsigned int) qrySTPtr->seqLenSL;
   } /*If: need more memory for sequence*/

   samSTPtr->readLenUI =
      cpStr_ulCp(
         samSTPtr->seqStr,
         qrySTPtr->seqStr
      );

   indexToSeq_alnSet(samSTPtr->seqStr);
     /*sequence is in index format, but I need a human
     `   readable format
     */

   /*****************************************************\
   * Fun34 Sec05 Sub04:
   *   - copy q-score
   \*****************************************************/

   if(
         qrySTPtr->qStr
      && qrySTPtr->qStr[0] != '\0'
   ){ /*If: have q-score entry*/

      if(samSTPtr->qSizeUI < qrySTPtr->seqLenSL)
      { /*If: need more memory*/
         free(samSTPtr->qStr);
         samSTPtr->qStr = 0;

         samSTPtr->qStr =
            malloc(
                 (qrySTPtr->seqLenSL + 9)
               * sizeof(unsigned char)
            );

         if(! samSTPtr)
            goto memErr_fun34_sec06;

         samSTPtr->seqSizeUI =
           (unsigned int) qrySTPtr->seqLenSL;
      } /*If: need more memory*/

      cpQEntry_samEntry(
         samSTPtr,
         qrySTPtr->qStr,
         0 /*histgram and q-score values already blanked*/
      ); /*copy q-score entries*/
   } /*If: have q-score entry*/

   else
   { /*Else: mark no q-score entry*/
      samSTPtr->qStr[0] = '*';
      samSTPtr->qStr[1] = '\0';
   } /*Else: mark no q-score entry*/

   /*****************************************************\
   * Fun34 Sec05 Sub05:
   *   - find mutation counts
   \*****************************************************/

   /*this would be faster if I did it as I aligned,
   `   however, this is eaiser. It should also add little
   `   overhead
   */

   for(
      siCig = 0;
      siCig < (signed int) samSTPtr->cigLenUI;
      ++siCig
   ){ /*Loop: get mutation counts*/

      switch(samSTPtr->cigTypeStr[siCig])
      { /*Switch: find mutation type*/
         case '=':
         /*Case: match*/
            samSTPtr->matchCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         /*Case: match*/

         case 'X':
         /*Case: mismatch*/
            samSTPtr->snpCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         /*Case: mismatch*/

         case 'I':
         /*Case: insertion*/
            samSTPtr->insCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         /*Case: insertion*/

         case 'D':
         /*Case: deletion*/
            samSTPtr->delCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         /*Case: deletion*/

         case 'S':
         /*Case: soft mask*/
            samSTPtr->maskCntUI +=
               samSTPtr->cigArySI[siCig];
            break;
         /*Case: soft mask*/
      } /*Switch: find mutation type*/

   }  /*Loop: get mutation counts*/

   samSTPtr->alnReadLenUI = samSTPtr->matchCntUI;
   samSTPtr->alnReadLenUI += samSTPtr->snpCntUI;
   samSTPtr->alnReadLenUI += samSTPtr->delCntUI;

   samSTPtr->refEndUI = samSTPtr->refStartUI;
   samSTPtr->refEndUI += samSTPtr->alnReadLenUI;

   if(samSTPtr->cigTypeStr[samSTPtr->cigLenUI - 1] != 'S')
      ;
   else if(samSTPtr->cigArySI[samSTPtr->cigLenUI -1] == 1)
   { /*If: need to make sure did not overshoot*/
      siCig = samSTPtr->matchCntUI;
      siCig += samSTPtr->snpCntUI;
      siCig += samSTPtr->insCntUI;
      siCig += samSTPtr->maskCntUI;

      if(siCig > (signed int) samSTPtr->readLenUI)
      { /*If: overshot*/
         --samSTPtr->maskCntUI;
         --samSTPtr->cigLenUI;
         samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = 0;
         samSTPtr->cigArySI[samSTPtr->cigLenUI] = 0;
      } /*If: overshot*/
   } /*If: need to make sure did not overshoot*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun34 Sec06:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun34_sec06;

   memErr_fun34_sec06:;
      errSC = def_memErr_mapRead;
      goto ret_fun34_sec06;
   ret_fun34_sec06:;
      qrySTPtr->offsetSL = 0;
      qrySTPtr->endAlnSL = qrySTPtr->seqLenSL - 1;

      refSTPtr->offsetSL = 0;
      refSTPtr->endAlnSL = refSTPtr->seqLenSL - 1;

      if(cigHeapArySI)
         free(cigHeapArySI);
      cigHeapArySI = 0;

      if(cigHeapStr)
         free(cigHeapStr);
      cigHeapStr = 0;

      return errSC;
} /*mergeToSam_mapRead*/

/*-------------------------------------------------------\
| Fun35: scoreSubAln_mapRead
|   - scores the aligment in samEntry and finds best
|     sub-aligment
| Input:
|   - samSTPtr:
|     o samEntry struct pointer with alignment to score
|     o do not trim
|   - subBl:
|     o modifiy alignment to best sub-alignment
|       * 1: look for best sub-alignment
|       * 0: do not look for sub-alignments
|   - qrySTPtr:
|     o seqST struct pointer with query sequence used in
|       alignment
|     o seqStr must be as lookup index (seqToIndex_alnSet)
|   - refSTPtr:
|     o seqST struct pointer with reference sequence used
|       in alignment
|     o seqStr must be as lookup index (seqToIndex_alnSet)
|   - alnSetSTPTr:
|     o alnSet struct pointer with scoring settings
|   - maxScoreSLPtr:
|     o signed long pointer to hold maximum score possible
|       for the aligned part of the query
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o extraStr in samSTPtr to have score "AS:i:score"
|       entry
|     o cigTypeStr and cigArySI in samSTPtr to be setup
|       for best sub-alignment (if subBl is 1)
|     o cigLenUI in samSTPtr to be new cigar length if
|       looking for best sub-alignment
|     o maxScoreSLPtr to have maximum possible score for
|       returned alignment (exlcudes unaligned bases)
|     o errSCPtr:
|       * 0 for no errors
|       * def_noScore_mapRead if no sub-score greater
|         than 0 (unlikley)
|       * def_memErr_mapRead for memory errors
|   - Returns:
|     o score for best sub-alignment or if subSL is 0, for
|       the entire samEntry
\-------------------------------------------------------*/
signed long
scoreSub_mapRead(
   struct samEntry *samSTPtr, /*alignment to score*/
   signed char subBl,         /*1:look for sub-alingment*/
   struct seqST *qrySTPtr,    /*query sequence*/
   struct seqST *refSTPtr,    /*reference sequence*/
   struct alnSet *alnSetSTPtr,/*scoring settings*/
   signed long *maxScoreSLPtr,/*gets max possible score*/
   signed char *errSCPtr      /*errors*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun35 TOC:
   '   - scores the aligment in samEntry and finds best
   '     sub-aligment
   '   o fun35 sec01:
   '     - variable declartions
   '   o fun35 sec02:
   '     - allocate memory for assembler score (if need)
   '   o fun35 sec03:
   '     - find score and best sub-score
   '   o fun35 sec04:
   '     - convert cigar to sub-alingment cigar; if subaln
   '   o fun35 sec05:
   '     - get extra entries and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun35 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *refStr = refSTPtr->seqStr;
   signed char *qryStr = qrySTPtr->seqStr;

   signed char *swapSCPtr = 0;
   unsigned int uiCig = 0;
   unsigned int uiPos = 0;
   signed long scoreSL = 0;

   signed long subScoreSL = 0;
   unsigned int subStartUI = 0;
   signed long maxSubScoreSL = 0;

   signed long bestScoreSL = 0;
   unsigned int bestStartUI = 0;
   unsigned int bestEndUI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun35 Sec02:
   ^   - allocate memory for assembler score (if need)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
      samSTPtr->extraLenUI >= samSTPtr->extraSizeUI -32
   ){ /*If: need more memory*/
      /*32 is more then enough for two numeric entries*/

      swapSCPtr =
         realloc(
            samSTPtr->extraStr,
              (samSTPtr->extraLenUI + 41)
            * sizeof(signed char)
         );
      if(! swapSCPtr)
         goto memErr_fun35_sec05;
      samSTPtr->extraStr = swapSCPtr;

      samSTPtr->extraSizeUI += 32;
   }  /*If: need more memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun35 Sec03:
   ^   - find score and best sub-score
   ^   o fun35 sec03 sub01:
   ^     - get to first aligned reference base
   ^   o fun35 sec03 sub02:
   ^     - score matches/snps + start switch + start loop
   ^   o fun35 sec03 sub03:
   ^     - score insertions
   ^   o fun35 sec03 sub04:
   ^     - score deletions
   ^   o fun35 sec03 sub05:
   ^     - masking, move to next query entry
   ^   o fun35 sec03 sub06:
   ^     - add final sub-score to total score
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun35 Sec03 Sub01:
   *   - get to first aligned reference base
   \*****************************************************/

   refStr += samSTPtr->refStartUI;
   *maxScoreSLPtr = 0;

   /*****************************************************\
   * Fun35 Sec03 Sub02:
   *   - score matches/snps + start switch + start loop
   \*****************************************************/

   for(
      uiCig = 0;
      uiCig < samSTPtr->cigLenUI;
      ++uiCig
   ){ /*Loop: score cigar*/

      switch(samSTPtr->cigTypeStr[uiCig])
      { /*Switch: find cigar type*/

         case '=':
         case 'X':
         case 'M':
         /*Case: match or snp*/
            for(
               uiPos = 0;
               uiPos <
                 (unsigned int) samSTPtr->cigArySI[uiCig];
               ++uiPos
            ){ /*Loop: find score*/
               subScoreSL +=
                  getScore_alnSet(
                     *qryStr,
                     *refStr,
                     alnSetSTPtr
                  );

               maxSubScoreSL +=
                  getScore_alnSet(
                     *qryStr,
                     *qryStr,
                     alnSetSTPtr
                  );

               ++refStr;
               ++qryStr;
            }  /*Loop: find score*/

            if(subScoreSL <= 0)
            { /*If: end of sub-alignment*/
               scoreSL += subScoreSL;
               *maxScoreSLPtr += maxSubScoreSL;

               subScoreSL = 0;
               maxSubScoreSL = 0;

               subStartUI = uiCig;
            } /*If: end of sub-alignment*/

            else if(subScoreSL > bestScoreSL)
            { /*Else If: new best score*/
               bestScoreSL = subScoreSL;
               bestStartUI = subStartUI;
               bestEndUI = uiCig;
            } /*Else If: new best score*/

            break;
         /*Case: match or snp*/

         /***********************************************\
         * Fun35 Sec03 Sub03:
         *   - score insertions
         \************************************************/

         case 'I':
         /*Case: insertion*/

            /*get score and subscore*/
            subScoreSL += alnSetSTPtr->gapSS;

            subScoreSL +=
               (
                    (samSTPtr->cigArySI[uiCig] - 1)
                  * alnSetSTPtr->extendSS
               );

            /*get the maximum possible score*/
            for(
               uiPos = 0;
               uiPos <
                 (unsigned int) samSTPtr->cigArySI[uiCig];
               ++uiPos
            ){ /*Loop: find maximum possible score*/
               maxSubScoreSL +=
                  getScore_alnSet(
                     *qryStr,
                     *qryStr,
                     alnSetSTPtr
                  );

               ++qryStr;
            }  /*Loop: find maximum possible score*/


            /*check if end of sub-alignment*/
            if(subScoreSL <= 0)
            { /*If: end of sub-alignment*/
               scoreSL += subScoreSL;
               *maxScoreSLPtr += maxSubScoreSL;

               subScoreSL = 0;
               maxSubScoreSL = 0;

               subStartUI = uiCig;
            } /*If: end of sub-alignment*/


            break;
         /*Case: insertion*/

         /***********************************************\
         * Fun35 Sec03 Sub04:
         *   - score deletions
         \************************************************/

         case 'D':
         /*Case: deletion*/
            refStr += samSTPtr->cigArySI[uiCig];
            subScoreSL += alnSetSTPtr->gapSS;

            subScoreSL +=
               (
                    (samSTPtr->cigArySI[uiCig] - 1)
                  * alnSetSTPtr->extendSS
               );

            if(subScoreSL <= 0)
            { /*If: end of sub-alignment*/
               scoreSL += subScoreSL;
               *maxScoreSLPtr += maxSubScoreSL;

               subScoreSL = 0;
               maxSubScoreSL = 0;

               subStartUI = uiCig;
            } /*If: end of sub-alignment*/

            break;
         /*Case: deletion*/

         /***********************************************\
         * Fun35 Sec03 Sub05:
         *   - masking, move to next query entry
         \************************************************/

         case 'S':
            qryStr += samSTPtr->cigArySI[uiCig];
            break;
         /*else is hard maksed*/
      } /*Switch: find cigar type*/

   } /*Loop: score cigar*/

   /*****************************************************\
   * Fun35 Sec03 Sub06:
   *   - add final sub-score to total score
   \*****************************************************/

   scoreSL += subScoreSL;
   *maxScoreSLPtr += maxSubScoreSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun35 Sec04:
   ^   - convert cigar to sub-alingment cigar (if subaln)
   ^   o fun35 sec04 sub01:
   ^     - check if have subalinment
   ^   o fun35 sec04 sub02:
   ^     - adjusting cigar end for sub-alignment
   ^   o fun35 sec04 sub03:
   ^     - handel adjusting start for sub-alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun35 Sec04 Sub01:
   *   - check if have sub-alignment
   \*****************************************************/

   if( (bestEndUI + 1) >= samSTPtr->cigLenUI )
      ; /*no softmasking at end or at end*/
   else if(samSTPtr->cigTypeStr[bestEndUI + 1] == 'S')
      ++bestEndUI;
      /*softmasking at aligment, so keeping end of aln*/

   if(! bestStartUI)
      ; /*subalignment starts at first base*/
   else if(bestStartUI > 1)
      ; /*need to mask to get best sub-alignment*/
   else if(samSTPtr->cigTypeStr[0] == 'S')
      bestStartUI = 0; /*no change to alignment*/

   if(! subBl)
      goto ret_fun35_sec05;
   else if(bestStartUI == bestEndUI)
      goto noScore_fun35_sec05;
   else if(bestStartUI >= samSTPtr->cigLenUI)
      goto noScore_fun35_sec05;

   /*****************************************************\
   * Fun35 Sec04 Sub02:
   *   - adjusting cigar end for sub-alignment
   \*****************************************************/

   if(bestEndUI + 1 < samSTPtr->cigLenUI)
   { /*If: removing ending cigar entries*/
      uiCig = bestEndUI + 1;
      subScoreSL = 0; /*using to count masked bases*/
     
      while(uiCig < samSTPtr->cigLenUI)
      { /*Loop: remove tailing entries*/

         switch(samSTPtr->cigTypeStr[uiCig])
         { /*Switch: find cigar type*/

            case 'M':
            case '=':
            /*Case: match or match/snp*/
               subScoreSL += samSTPtr->cigArySI[uiCig];

               samSTPtr->matchCntUI -=
                 samSTPtr->cigArySI[uiCig];
               samSTPtr->refEndUI -=
                  samSTPtr->cigArySI[uiCig];
               samSTPtr->alnReadLenUI -=
                  samSTPtr->cigArySI[uiCig];
               break;
            /*Case: match or match/snp*/

            case 'X':
            /*Case: snp*/
               subScoreSL += samSTPtr->cigArySI[uiCig];

               samSTPtr->snpCntUI -=
                 samSTPtr->cigArySI[uiCig];
               samSTPtr->refEndUI -=
                  samSTPtr->cigArySI[uiCig];
               samSTPtr->alnReadLenUI -=
                  samSTPtr->cigArySI[uiCig];
               break;
            /*Case: snp*/

            case 'D':
            /*Case: deletion*/
               samSTPtr->delCntUI -=
                 samSTPtr->cigArySI[uiCig];
               samSTPtr->refEndUI -=
                  samSTPtr->cigArySI[uiCig];
               samSTPtr->alnReadLenUI -=
                  samSTPtr->cigArySI[uiCig];
               break;
            /*Case: deletion*/

            case 'I':
            /*Case: insertion*/
               subScoreSL += samSTPtr->cigArySI[uiCig];

               samSTPtr->insCntUI -=
                 samSTPtr->cigArySI[uiCig];
               break;
            /*Case: insertion*/

            case 'S':
               subScoreSL += samSTPtr->cigArySI[uiCig];
               break;
         } /*Switch: find cigar type*/

         ++uiCig;
      } /*Loop: remove tailing entries*/

      uiCig = bestEndUI + 1;
      samSTPtr->cigTypeStr[uiCig] = 'S';
      samSTPtr->cigArySI[uiCig] = (signed int) subScoreSL;
      samSTPtr->cigLenUI = uiCig + 1;
      samSTPtr->cigTypeStr[samSTPtr->cigLenUI] = '\0';
   } /*If: removing ending cigar entries*/


   /*****************************************************\
   * Fun35 Sec04 Sub03:
   *   - handel adjusting start for sub-alignment
   *   o fun35 sec04 sub03 cat01:
   *     - start removal; convert first cigar entry to S
   *   o fun35 sec04 sub03 cat02:
   *     - start removal; remove targed cigar entries
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun35 Sec04 Sub03 Cat01:
   +   - start removal; convert first cigar entry to S
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(bestStartUI)
   { /*If: need to shift cigar*/
      if(
            samSTPtr->cigTypeStr[0] != 'S'
         && samSTPtr->cigTypeStr[0] != 'H'
      ){ /*If: start not masked*/
         samSTPtr->maskCntUI += samSTPtr->cigArySI[0];

         switch(samSTPtr->cigTypeStr[0])
         { /*Switch: find cigar type*/

            case 'M':
            case '=':
            /*Case: match or match/snp*/
               samSTPtr->matchCntUI -=
                 samSTPtr->cigArySI[0];
               samSTPtr->refStartUI +=
                 samSTPtr->cigArySI[0];
               samSTPtr->alnReadLenUI -=
                 samSTPtr->cigArySI[0];
               break;
            /*Case: match or match/snp*/

            case 'X':
            /*Case: snp*/
               samSTPtr->snpCntUI -=
                 samSTPtr->cigArySI[0];
               samSTPtr->refStartUI +=
                 samSTPtr->cigArySI[0];
               samSTPtr->alnReadLenUI -=
                 samSTPtr->cigArySI[0];
               break;
            /*Case: snp*/

            case 'D':
            /*Case: deletion*/
               samSTPtr->delCntUI -=
                 samSTPtr->cigArySI[0];

               samSTPtr->refStartUI +=
                 samSTPtr->cigArySI[0];
               samSTPtr->alnReadLenUI -=
                 samSTPtr->cigArySI[0];

               samSTPtr->cigArySI[0] = 0;
               break;
            /*Case: deletion*/

            case 'I':
            /*Case: insertion*/
               samSTPtr->insCntUI -=
                 samSTPtr->cigArySI[0];
               break;
            /*Case: insertion*/
         } /*Switch: find cigar type*/

      } /*If: start not masked*/

      samSTPtr->cigTypeStr[0] = 'S';

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun35 Sec04 Sub03 Cat02:
      +   - start removal; remove targed cigar entries
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      uiCig = 1;
      uiPos = 1;
     
      while(uiCig < samSTPtr->cigLenUI)
      { /*Loop: shift cigar*/

         if(uiCig <= bestStartUI)
         { /*If: masking*/

            if(samSTPtr->cigTypeStr[uiCig] != 'D')
            { /*If: not a deletion*/
               samSTPtr->cigArySI[0] +=
                  samSTPtr->cigArySI[uiCig];
               samSTPtr->maskCntUI +=
                  samSTPtr->cigArySI[uiCig];
            } /*If: not a deletion*/

            if(samSTPtr->cigTypeStr[uiCig] != 'I')
            { /*If: not an insertion*/
               samSTPtr->refStartUI +=
                 samSTPtr->cigArySI[uiCig];
               samSTPtr->alnReadLenUI -=
                 samSTPtr->cigArySI[uiCig];
            } /*If: not an insertion*/

            switch(samSTPtr->cigTypeStr[uiCig])
            { /*Switch: find cigar type*/

               case 'M':
               case '=':
                  samSTPtr->matchCntUI -=
                    samSTPtr->cigArySI[uiCig];
                  break;

               case 'X':
                  samSTPtr->snpCntUI -=
                    samSTPtr->cigArySI[uiCig];
                  break;

               case 'D':
                  samSTPtr->delCntUI -=
                    samSTPtr->cigArySI[uiCig];
                  break;

               case 'I':
                  samSTPtr->insCntUI -=
                    samSTPtr->cigArySI[uiCig];
                  break;
            } /*Switch: find cigar type*/
         } /*If: masking*/

         samSTPtr->cigTypeStr[uiPos] =
            samSTPtr->cigTypeStr[uiCig];
         samSTPtr->cigArySI[uiPos] =
            samSTPtr->cigArySI[uiCig];

         if(uiCig > bestStartUI)
            ++uiPos;

         ++uiCig;
      }  /*Loop: shift cigar*/

      samSTPtr->cigLenUI = uiPos;
      samSTPtr->cigTypeStr[uiPos] = '\0';
      bestEndUI = samSTPtr->cigLenUI;
   } /*If: need to shift cigar*/

   /*****************************************************\
   * Fun35 Sec04 Sub04:
   *   - adjusting alignment lengths and coordiantes
   \*****************************************************/

   samSTPtr->alnReadLenUI = samSTPtr->matchCntUI;
   samSTPtr->alnReadLenUI += samSTPtr->snpCntUI;
   samSTPtr->alnReadLenUI += samSTPtr->delCntUI;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun35 Sec05:
   ^   - get extra entries and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   ret_fun35_sec05:;
      *errSCPtr = 0;

      samSTPtr->extraStr[samSTPtr->extraLenUI++] = 'N';
      samSTPtr->extraStr[samSTPtr->extraLenUI++] = 'M';
      samSTPtr->extraStr[samSTPtr->extraLenUI++] = ':';
      samSTPtr->extraStr[samSTPtr->extraLenUI++] = 'i';
      samSTPtr->extraStr[samSTPtr->extraLenUI++] = ':';

      samSTPtr->extraLenUI +=
         numToStr(
            &samSTPtr->extraStr[samSTPtr->extraLenUI],
            (unsigned long)
            (
                samSTPtr->snpCntUI
              + samSTPtr->delCntUI
              + samSTPtr->insCntUI
            )
         );

      samSTPtr->extraStr[samSTPtr->extraLenUI++] = '\t';

      samSTPtr->extraStr[samSTPtr->extraLenUI++] = 'A';
      samSTPtr->extraStr[samSTPtr->extraLenUI++] = 'S';
      samSTPtr->extraStr[samSTPtr->extraLenUI++] = ':';
      samSTPtr->extraStr[samSTPtr->extraLenUI++] = 'i';
      samSTPtr->extraStr[samSTPtr->extraLenUI++] = ':';

      bestScoreSL = scoreSL; /*so can  update AS:i score*/

      if(scoreSL < 0)
      { /*If: had negative score*/
         samSTPtr->extraStr[samSTPtr->extraLenUI++] = '-';
         bestScoreSL *= -1;
      } /*If: had negative score*/
      
      samSTPtr->extraLenUI +=
         numToStr(
            &samSTPtr->extraStr[samSTPtr->extraLenUI],
            (unsigned long) bestScoreSL
         );

      return scoreSL;

   noScore_fun35_sec05:;
      *errSCPtr = def_noScore_mapRead;
      return scoreSL;

   memErr_fun35_sec05:;
      *errSCPtr = def_memErr_mapRead;
      return scoreSL;
} /*scoreSubAln_mapRead*/

/*-------------------------------------------------------\
| Fun36: chainToAln_mapRead
|   - converts chain in chains_mapRead struct to alignment
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to add alignment to
|   - chainsSTPtr:
|     o chains_mapRead struct pointer with chains to use
|       in making the alignment
|   - indexSI:
|     o starting index of merged chain to align
|   - qrySTPtr:
|     o seqST struct with query sequence
|   - refSTPtr:
|     o seqST struct pointer with reference sequence
|   - setSTPtr:
|     o set_mapRead struct pointer with settings to use
|   - lenKmerUC:
|     o length of one kmer in reference/query
|   - chainLenSI:
|     o minimum length of a chain
|     o used to align ends before chains, alignment size
|       is 1.5 chains (for chainLenSI = 20, is 30 bases)
|       * 1.5 is off the top of my head
|   - matrixSTPtr:
|     o directional matrix to use in alignment step
|   - errSCPtr:
|     o signed char pointer to hold errors
| Output:
|   - Modifies:
|     o samSTPtr to have alignment (if could align)
|     o uses matrixSTPtr for aligment step, is modified,
|       and resized, but nothing of relavence
|     o sets errSCPtr
|       * 0 for no errors
|       * def_memErr_mapRead for memory errors
|       * def_noChains_mapRead if merged chains at index
|         did not have enough coverage
|       * def_noAln_mapRead if aligment was under
|         minimum percent query length or score
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|     o def_noChains_mapRead if merged chains at index did
|       not have enough coverage
|     o def_noAln_mapRead if aligment was under
|       minimum percent query length or score
\-------------------------------------------------------*/
signed long
chainToAln_mapRead(
   struct samEntry *samSTPtr,    /*gets alignment*/
   struct chains_mapRead *chainsSTPtr, /*chains to align*/
   signed int indexSI,           /*chain index to align*/
   struct seqST *qrySTPtr,       /*query sequence*/
   struct seqST *refSTPtr,       /*reference sequence*/
   struct set_mapRead *setSTPtr, /*alignment settings*/
   unsigned char lenKmerUC,      /*lenght of 1 kemr*/
   signed int chainLenSI,        /*min length of a chain*/
   struct dirMatrix *matrixSTPtr,/*for alignment*/
   signed char *errSCPtr         /*gets errors*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun36 TOC:
   '   - converts chain in chains_mapRead struct to
   '     an alignment
   '   o fun36 sec01:
   '     - variable declarations
   '   o fun36 sec02:
   '     - align chains
   '   o fun36 sec03:
   '     - score alignment & check if passes min settings
   '   o fun36 sec04:
   '     - set error type and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun36 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   float percF = 0;
   signed long scoreSL = 0;
   signed long maxScoreSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun36 Sec02:
   ^   - align chains
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_samEntry(samSTPtr);

   /*make sure this alignment could be kept*/
   percF = (float) chainsSTPtr->scoreArySI[indexSI];
   percF += (float) lenKmerUC;
   percF -= 2;
      /*-1 for index 1 to 0 coversion, another -1 (-2)
      `  to account for one base always kept
      */
   percF /= (float) qrySTPtr->seqLenSL;

   if(percF < setSTPtr->chainMinNtF)
      goto noChains_fun36_sec04;

   if(
      mergeToSam_mapRead(
         chainsSTPtr,
         indexSI,
         setSTPtr->alnEndsBl,
         qrySTPtr,
         refSTPtr,
         lenKmerUC,
         chainLenSI,
         setSTPtr->alnSetST,
         matrixSTPtr,
         samSTPtr
      ) /*make merged chains into an alignment*/
   ) goto memErr_fun36_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun36 Sec03:
   ^   - score alignment and check if passes min settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   scoreSL =
      scoreSub_mapRead(
         samSTPtr,
         setSTPtr->subBl,  /*if doing sub-alignment*/
         qrySTPtr,
         refSTPtr,
         setSTPtr->alnSetST,
         &maxScoreSL,     /*gets maximum score possible*/
         errSCPtr
      );


   /*check if had the minimum percent score*/
   percF = (float) scoreSL;
   percF /= (float) maxScoreSL;

   if(percF < setSTPtr->minScoreF)
      goto noAln_fun36_sec04;


   /*check if had the minimum percent length*/
   percF = (float) samSTPtr->matchCntUI;
   percF += (float) samSTPtr->snpCntUI;
   percF += (float) samSTPtr->insCntUI;
   percF /= (float) qrySTPtr->seqLenSL;

   if(percF < setSTPtr->minPercLenF)
      goto noAln_fun36_sec04;


   /*find percent matches for aligned read*/
   percF = (float) samSTPtr->matchCntUI;

   maxScoreSL = (signed long) samSTPtr->matchCntUI;
   maxScoreSL += (signed long) samSTPtr->snpCntUI;
   maxScoreSL += (signed long) samSTPtr->insCntUI;

   percF /= (float) maxScoreSL;

   if(percF < setSTPtr->minMatchF)
      goto noAln_fun36_sec04;
      /*not enough matches*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun36 Sec04:
   ^   - set error type and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *errSCPtr = 0;
   goto ret_fun36_sec04;

   memErr_fun36_sec04:;
      *errSCPtr = def_memErr_mapRead;
      goto ret_fun36_sec04;

   noChains_fun36_sec04:;
      *errSCPtr = def_noChains_mapRead;
      goto ret_fun36_sec04;

   noAln_fun36_sec04:;
      *errSCPtr = def_noAln_mapRead;
      goto ret_fun36_sec04;

   ret_fun36_sec04:;
      return scoreSL;
} /*chainToAln_mapRead*/

/*-------------------------------------------------------\
| Fun37: align_mapRead
|   - maps read to reference
| Input:
|   - qrySTPtr:
|     o seqST struct pointer with query sequence and id
|   - refSTPtr:
|     o ref_mapRead struct pointer with reference
|       sequence, settings, and kmer arrays for alignment
|     o setup with init_ref_mapRead, setup_ref_mapRead,
|       and addRef_ref_mapRead
|   - mapIndexSIPtr:
|     o signed int pointer to hold index of best chain
|       (used to make mapping)
|   - alnSTPtr:
|     o aln_mapRead struct pionter with allocated memory to
|       use/allocate/resize in mapping steps
|   - setSTPtr:
|     o set_mapRead structure with non-reference specific
|   - samSTPtr:
|     o pointer to samEntry struct to hold best alignment
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o samSTPtr to have highest scoring alignment
|     o seqStr in qrySTPtr to be lookup index's, you can
|       undo this with indexToSeq_alnSet(qrySTPtr->seqStr)
|     o mapIndexSLPtr to have index of chain used to find
|       the best alignment
|     o kmerArySI in alnSTPtr to have query kmer sequence
|     o lenSI in alnSTPtr to be length of kmerArySI
|     o sizeSI in alnSTPtr if kmerArySI is resized
|     o chainsSTPtr alnSTPtr to have found chains
|     o errSCPtr to have error message
|       * 0 for no errors
|       * def_memErr_mapRead for memory errors
|       * def_noChains_mapRead if chould not find any
|         chains
|       * def_noAln_mapRead if aligment was under
|         minimum percent query length or score
|     o matrixSTPtr in alnSTPtr to nothing usefull. Here
|       to avoid memory allocations in alignment step
|   - Returns:
|     o score for primary alignment
\-------------------------------------------------------*/
signed long
align_mapRead(
   struct seqST *qrySTPtr,      /*has query sequence*/
   struct ref_mapRead *refSTPtr,/*ref kmers/settings*/
   signed int *mapIndexSIPtr,   /*index of best chain*/
   struct samEntry *samSTPtr,   /*gets alignments*/
   struct aln_mapRead *alnSTPtr,/*memory for query*/
   struct set_mapRead *setSTPtr,/*alignment settings*/
   signed char *errSCPtr        /*gets errors*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun37 TOC:
   '   - maps a read to a reference
   '   o fun31 sec01:
   '     - variable declarations
   '   o fun37 sec02:
   '     - deal with to short alignments (waterman)
   '   o fun37 sec03:
   '     - setup and get query kmer sequence
   '   o fun37 sec04:
   '     - get chains, merge chains, and check coverage
   '   o fun37 sec05:
   '     - get alignments
   '   o fun37 sec06:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun37 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed long scoreSL = 0;
   signed long bestScoreSL = def_minSL_mapRead;
   unsigned int ignoreUI = 0;

   float percF = 0;

   signed int siIndex = 0;
   struct hit_mapRead *hitHeapST = 0;
   struct hit_mapRead *tailSTPtr = 0;
   struct hit_mapRead *iterSTPtr = 0;
   struct hit_mapRead *swapSTPtr = 0;

   struct samEntry samStackST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun37 Sec02:
   ^   - deal with to short alignments (waterman)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_samEntry(&samStackST);
   if( setup_samEntry(&samStackST) )
      goto memErr_fun37_sec06;

   seqToIndex_alnSet(qrySTPtr->seqStr);

   if(qrySTPtr->seqLenSL < (refSTPtr->lenKmerUC << 1))
   { /*If: mapping a read that has less then two chains*/
      qrySTPtr->offsetSL = 0;
      qrySTPtr->endAlnSL = qrySTPtr->seqLenSL - 1;

      refSTPtr->seqSTPtr->offsetSL = 0;
      refSTPtr->seqSTPtr->endAlnSL =
         refSTPtr->seqSTPtr->seqLenSL - 1;

      bestScoreSL =
         water(
            qrySTPtr,
            refSTPtr->seqSTPtr,
            alnSTPtr->matrixSTPtr,
            setSTPtr->alnSetST
         );

      if(alnSTPtr->matrixSTPtr->errSC == def_memErr_water)
         goto memErr_fun37_sec06;
      if(
         getAln_dirMatrix(
            alnSTPtr->matrixSTPtr,
            0, /*0 = use index in matrixSTPtr*/
            0, /*forward alignment*/
            qrySTPtr,
            refSTPtr->seqSTPtr,
            samSTPtr,
            &ignoreUI,
            setSTPtr->alnSetST
         )
      ) goto memErr_fun37_sec06;

      scoreSL =
         maxScore_alnSet(
            qrySTPtr->seqStr,
            0,
            qrySTPtr->seqLenSL,
            setSTPtr->alnSetST
        );

      /*check if had the minimum percent score*/
      percF = (float) bestScoreSL;
      percF /= (float) scoreSL;

      if(percF < setSTPtr->minScoreF)
         goto noAln_fun37_sec06;

      /*check if had the minimum percent length*/
      percF = (float) samSTPtr->matchCntUI;
      percF += (float) samSTPtr->snpCntUI;
      percF += (float) samSTPtr->insCntUI;
      percF /= (float) qrySTPtr->seqLenSL;

      if(percF < setSTPtr->minPercLenF)
         goto noAln_fun37_sec06;


      /*find percent matches for aligned read*/
      percF = (float) samSTPtr->matchCntUI;

      scoreSL = (signed long) samSTPtr->matchCntUI;
      scoreSL += (signed long) samSTPtr->snpCntUI;
      scoreSL += (signed long) samSTPtr->insCntUI;

      percF /= (float) scoreSL;

      if(percF < setSTPtr->minMatchF)
         goto noAln_fun37_sec06;
         /*not enough matches*/

      goto done_fun37_sec06;
   } /*If: mapping a read that has less then two chains*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun37 Sec03:
   ^   - setup and get query kmer sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   alnSTPtr->lenSI =
      seqToKmer_kmerFun(
         qrySTPtr->seqStr,
         qrySTPtr->seqLenSL,
         &alnSTPtr->kmerArySI,
         &alnSTPtr->sizeSI,
         refSTPtr->lenKmerUC
      );

   if(! alnSTPtr->lenSI)
      goto memErr_fun37_sec06;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun37 Sec04:
   ^   - get chains, merge chains, and check coverage
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
      kmerChain_mapRead(
         alnSTPtr->chainsSTPtr,
         refSTPtr->minChainLenSI,
         refSTPtr->lenKmerUC,
         alnSTPtr->kmerArySI,
         alnSTPtr->lenSI,
         refSTPtr->kmerArySI,
         refSTPtr->indexArySI,
         refSTPtr->lenSI,
         0
      )
   ) goto memErr_fun37_sec06;

   *mapIndexSIPtr =
       mergeChains_mapRead(
          alnSTPtr->chainsSTPtr,
          refSTPtr->lenKmerUC,
          setSTPtr->gapSI,    /*score for having gaps*/
          setSTPtr->matchSI,  /*score for match*/
          (signed int)
             (
                 qrySTPtr->seqLenSL 
               * setSTPtr->maxLenPercF
             )
             /*skip merge if gap is 25% of query length*/
       ); /*merges chains to get an alignment guide*/

   if(*mapIndexSIPtr == def_noChains_mapRead)
      goto noChains_fun37_sec06;

   /*make sure this alignment could be kept*/
   percF =
      (float)
      alnSTPtr->chainsSTPtr->scoreArySI[*mapIndexSIPtr];
   percF += (float) refSTPtr->lenKmerUC;
   percF -= 2;
      /*-1 for index 1 to 0 coversion, another -1 (-2)
      `  to account for one base always kept
      */
   percF /= (float) qrySTPtr->seqLenSL;

   if(percF < setSTPtr->chainMinNtF)
      goto noAln_fun37_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun37 Sec05:
   ^   - get alignments
   ^   o fun37 sec03 sub01:
   ^     - check if already aligned index + start loop
   ^   o fun37 sec03 sub02:
   ^     - add index's path to no-align list
   ^   o fun37 sec03 sub03:
   ^     - get alingment for index
   ^   o fun37 sec03 sub04:
   ^     - score alignment and check if keeping
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun37 Sec03 Sub01:
   *   - check if already aligned index + start loop
   \*****************************************************/

   for(
      siIndex = 0;
      siIndex < alnSTPtr->chainsSTPtr->lenChainSI;
      ++siIndex
   ){ /*Loop: find best alignment*/

      if(! hitHeapST)
         ;
      else if(siIndex == hitHeapST->indexSI)
      { /*Else If: index is part of larger aln*/
         hitHeapST->indexSI =
            alnSTPtr->chainsSTPtr->nextArySI[siIndex];

         if(! hitHeapST->indexSI)
         { /*If: no new index to add*/
            iterSTPtr = hitHeapST->nextST;
            free(hitHeapST);
            hitHeapST = iterSTPtr;
            iterSTPtr = 0;

            if(! hitHeapST)
               tailSTPtr = 0;
         } /*If: no new index to add*/

         else if(! hitHeapST->nextST)
            ; /*no other index's*/

         else
         { /*Else: need to position new index*/
            iterSTPtr = hitHeapST;

            while(iterSTPtr->nextST)
            { /*Loop: find insert position*/
               if(
                    hitHeapST->indexSI
                  < iterSTPtr->nextST->indexSI
               ) break;

               iterSTPtr = iterSTPtr->nextST;
            } /*Loop: find insert position*/

            if(! iterSTPtr)
            { /*If: new index goes to end*/
               swapSTPtr = hitHeapST->nextST;
               tailSTPtr->nextST = hitHeapST;
               hitHeapST->nextST = 0;
               hitHeapST = swapSTPtr;
            } /*If: new index goes to end*/

            else
            { /*Else: index is inserted*/
              swapSTPtr = hitHeapST; 
              hitHeapST = hitHeapST->nextST;
              swapSTPtr->nextST = iterSTPtr->nextST;
              iterSTPtr->nextST = swapSTPtr;
            } /*Else: index is inserted*/
         } /*Else: need to position new index*/

         swapSTPtr = 0;
         iterSTPtr = 0;
         continue;
      } /*Else If: index is part of larger aln*/

      /**************************************************\
      * Fun37 Sec03 Sub02:
      *   - add index's path to no-align list
      \**************************************************/

      swapSTPtr = malloc( sizeof(struct hit_mapRead) );
      if(! swapSTPtr)
         goto memErr_fun37_sec06;
      swapSTPtr->indexSI =
         alnSTPtr->chainsSTPtr->nextArySI[siIndex];
      swapSTPtr->nextST = 0;

      if(! hitHeapST)
      { /*If: first hit*/
         hitHeapST = swapSTPtr;
         tailSTPtr = hitHeapST;
      } /*If: first hit*/

      else if(swapSTPtr->indexSI < hitHeapST->indexSI)
      { /*If: new index is new head*/
         swapSTPtr->nextST = hitHeapST;
         hitHeapST = swapSTPtr;
      } /*If: new index is new head*/

      else
      { /*Else: need to position new index*/
         iterSTPtr = hitHeapST->nextST;

         if(! iterSTPtr)
         { /*If: no second index*/
            hitHeapST->nextST = swapSTPtr;
            tailSTPtr = swapSTPtr;
         } /*If: no second index*/

         else
         { /*Else: have more than one index*/
            while(iterSTPtr->nextST)
            { /*Loop: find insert position*/
               if(
                    swapSTPtr->indexSI
                  < iterSTPtr->nextST->indexSI
               ) break;

               iterSTPtr = iterSTPtr->nextST;
            } /*Loop: find insert position*/

            if(! iterSTPtr)
            { /*If: index goes at end*/
               tailSTPtr->nextST = swapSTPtr;
               tailSTPtr = tailSTPtr->nextST;
            } /*If: index goes at end*/

            else
            { /*Else: index is inserted*/
              swapSTPtr->nextST = iterSTPtr->nextST; 
              iterSTPtr->nextST = swapSTPtr;
            } /*Else: index is inserted*/
         } /*Else: have more than one index*/
      } /*Else: need to position new index*/

      swapSTPtr = 0;
      iterSTPtr = 0;

      /**************************************************\
      * Fun37 Sec03 Sub03:
      *   - get alingment for index
      \**************************************************/

      scoreSL = 
         chainToAln_mapRead(
            &samStackST,
            alnSTPtr->chainsSTPtr,
            siIndex,
            qrySTPtr,
            refSTPtr->seqSTPtr,
            setSTPtr,
            refSTPtr->lenKmerUC,
            refSTPtr->minChainLenSI,
            alnSTPtr->matrixSTPtr,
            errSCPtr
         );

      if(*errSCPtr == def_memErr_mapRead)
          goto memErr_fun37_sec06;
      else if(*errSCPtr)
         continue; /*no alingment*/

      if(scoreSL > bestScoreSL)
      { /*If: new best alignment*/
         bestScoreSL = scoreSL;

         swap_samEntry(
            &samStackST,
            samSTPtr
         );
      } /*If: new best alignment*/
   } /*Loop: find best alignment*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun37 Sec06:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   done_fun37_sec06:;
      if(bestScoreSL == def_minSL_mapRead)
         goto noAln_fun37_sec06;

      *errSCPtr = 0;
   goto ret_fun37_sec06;

   memErr_fun37_sec06:;
      *errSCPtr = def_memErr_mapRead;
      goto err_fun37_sec06;

   noAln_fun37_sec06:;
      *errSCPtr = def_noAln_mapRead;
      goto err_fun37_sec06;

   noChains_fun37_sec06:;
      *errSCPtr = def_noChains_mapRead;
      goto err_fun37_sec06;

   err_fun37_sec06:;
      blank_samEntry(samSTPtr);
      goto ret_fun37_sec06;

   ret_fun37_sec06:;
      freeStack_samEntry(&samStackST);

      while(hitHeapST)
      { /*Loop: free structures*/
         swapSTPtr = hitHeapST->nextST;
         free(hitHeapST);
         hitHeapST = swapSTPtr;
      } /*Loop: free structures*/

      swapSTPtr = 0;
      iterSTPtr = 0;
      hitHeapST = 0;
      tailSTPtr = 0;

      return bestScoreSL;
} /*align_mapRead*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconveint / not possible, this code is under the
:   MIT license
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
