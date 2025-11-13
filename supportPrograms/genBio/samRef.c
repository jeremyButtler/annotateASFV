/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' samRef SOF: Start Of File
'   - structure functions to get references lengths from
'     a sam file
'   - keeping this separate from samEntry avoids usage
'     of strAry in other programs. This can save 4kb in
'     compile size
'   o header:
'     - Included libraries
'   o .h st01: refs_samRef
'     - holds reference length & name from sam file header
'   o fun01: blank_refs_samRef
'     - blanks a refs_samRef struct
'   o fun02: init_refs_samRef
'     - initializes a refs_samRef struct
'   o fun03: freeStack_refs_samRef
'     - frees variables in a refs_samRef struct
'   o fun04: freeHeap_refs_samRef
'     - frees a refs_samRef struct
'   o fun05: setup_refs_samRef
'     - allocates memory for a refs_samRef struct
'   o fun06: realloc_refs_samRef
'     - reallocates memory for a refs_samRef struct
'   o fun07: getRefLen_samRef
'     - gets reference ids & length from a sam file header
'   o fun08: findRef_refs_samRef
'     - finds a reference id in a refs_samRef struct
'   o fun09: addRef_samRef
'     - adds reference information to array in refStack
'   o fun10: buildRefMergeIndex_samRef
'     - looks for refseq ids that might be contigs from
'       the sam assembly and builds an index for each
'       unique id
'   o fun11: pSamHeader_samRef
'     - prints samEntry header for a reference or set of
'       references
'   o fun12: pSamToRef_samRef
'     - print samEntry to its reference bin
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

#include <stdio.h>

#include "samRef.h"
#include "samEntry.h"
#include "../genLib/strAry.h"
#include "../genLib/ulCp.h"
#include "../genLib/base10str.h"

/*These have no .c files*/
#include "../genLib/endLine.h" /*for OS line breaks*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden dependencies:
!   - .c  #include "../genLib/numToStr.h"
!   - .c  #include "../genLib/fileFun.h"
!   - .h  #include "ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: blank_refs_samRef
|   - blanks a refs_samRef struct
| Input:
|   - refsSTPtr:
|     o pointer to refs_samRef struct to blank
| Output:
|   - Modifies:
|     o refsSTPtr to be blanked
\-------------------------------------------------------*/
void
blank_refs_samRef(
   struct refs_samRef *refSTPtr
){
   if(refSTPtr->lenAryUI)
      refSTPtr->lenAryUI[0] = 0;

   if(refSTPtr->idAryStr)
      refSTPtr->idAryStr[0] = '\0';

   refSTPtr->numRefUI = 0;
} /*blank_refs_samRef*/

/*-------------------------------------------------------\
| Fun02: init_refs_samRef
|   - initializes a refs_samRef struct
| Input:
|   - refsSTPtr:
|     o pointer to refs_samRef struct to initialize
| Output:
|   - Modifies:
|     o refsSTPtr to be initialize (everything 0)
\-------------------------------------------------------*/
void
init_refs_samRef(
   struct refs_samRef *refSTPtr
){
   refSTPtr->lenAryUI = 0;
   refSTPtr->idAryStr = 0;
   refSTPtr->arySizeUI = 0;

   blank_refs_samRef(refSTPtr);
} /*init_refs_samRef*/

/*-------------------------------------------------------\
| Fun03: freeStack_refs_samRef
|   - frees variables in a refs_samRef struct
| Input:
|   - refsSTPtr:
|     o pointer to refs_samRef with variables to free
| Output:
|   - Frees:
|     o lenAryUI and idAryStr arrays in refsSTPtr
\-------------------------------------------------------*/
void
freeStack_refs_samRef(
   struct refs_samRef *refSTPtr
){
   if(! refSTPtr)
      return; /*no struct input*/

   if(refSTPtr->lenAryUI)
      free(refSTPtr->lenAryUI);

   if(refSTPtr->idAryStr)
      free(refSTPtr->idAryStr);

   init_refs_samRef(refSTPtr);
} /*freeStack_refs_samRef*/

/*-------------------------------------------------------\
| Fun04: freeHeap_refs_samRef
|   - frees a refs_samRef struct
| Input:
|   - refsSTPtr:
|     o pointer to refs_samRef struct to free
| Output:
|   - Frees:
|     o refSTPtr (you must set to 0)
\-------------------------------------------------------*/
void
freeHeap_refs_samRef(
   struct refs_samRef *refSTPtr
){
   if(refSTPtr)
   { /*If: have struct to free*/
      freeStack_refs_samRef(refSTPtr);
      free(refSTPtr);
   } /*If: have struct to free*/
} /*freeHeap_refs_samRef*/

/*-------------------------------------------------------\
| Fun05: setup_refs_samRef
|   - allocates memory for a refs_samRef struct
| Input:
|   - refsSTPtr:
|     o pointer to refs_samRef struct to get memory for
|   - numRefsUI:
|     o number of refs to start out with (0 goes to 16)
| Output:
|   - Mofidifies:
|     o lenAryUI and idAryUI in refSTPtr to have memory
|     o sizeAryUI to be numRefsUI
|     o calls blank_refs_samRef
|   - Returns:
|     o 0 for no errors
|     o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
setup_refs_samRef(
   struct refs_samRef *refSTPtr,
   unsigned int numRefsUI
){
   if(numRefsUI <= 0)
      numRefsUI = 16;

   if(refSTPtr->lenAryUI)
      free(refSTPtr->lenAryUI);
   refSTPtr->lenAryUI = 0;

   if(refSTPtr->idAryStr)
      free(refSTPtr->idAryStr);
   refSTPtr->idAryStr = 0;


   refSTPtr->lenAryUI =
      calloc(
         numRefsUI,
         sizeof(unsigned int)
      );

   if(! refSTPtr->lenAryUI)
      goto memErr_fun05;


   refSTPtr->idAryStr = mk_strAry(numRefsUI);

   if(! refSTPtr->idAryStr)
      goto memErr_fun05;

   refSTPtr->arySizeUI = numRefsUI;
   blank_refs_samRef(refSTPtr);

   return 0;

   memErr_fun05:;
   return def_memErr_samEntry;
} /*setup_refs_samRef*/

/*-------------------------------------------------------\
| Fun06: realloc_refs_samRef
|   - reallocates memory for a refs_samRef struct
| Input:
|   - refsSTPtr:
|     o pointer to refs_samRef struct to reallocate
|   - numRefsUI:
|     o new number of refs
| Output:
|   - Mofidifies:
|     o lenAryUI and idAryUI in refSTPtr to be resized
|     o sizeAryUI to be numRefsUI
|     o calls blank_refs_samRef
|   - Returns:
|     o 0 for no errors
|     o def_memErr_samEntry for memory errors
\-------------------------------------------------------*/
signed char
realloc_refs_samRef(
   struct refs_samRef *refSTPtr,
   unsigned int numRefsUI
){
   signed char *tmpSCPtr = 0;

   if(numRefsUI == refSTPtr->numRefUI)
      goto done_fun05; /*already correct size*/


   tmpSCPtr =
      (signed char *)
      realloc(
         refSTPtr->lenAryUI,
         numRefsUI * sizeof(unsigned int)
      );

   if(! tmpSCPtr)
      goto memErr_fun05;

   refSTPtr->lenAryUI = (unsigned int *) tmpSCPtr;


   tmpSCPtr =
      realloc_strAry(
         refSTPtr->idAryStr,
         numRefsUI
      );

   if(! tmpSCPtr)
      goto memErr_fun05;

   refSTPtr->idAryStr = tmpSCPtr;

   refSTPtr->arySizeUI = numRefsUI;

   done_fun05:;
      return 0;

   memErr_fun05:;
      return def_memErr_samEntry;
} /*realloc_refs_samRef*/

/*-------------------------------------------------------\
| Fun07: getRefLen_samRef
|   - gets reference ids and length from a sam file header
| Input:
|   - refSTPtr:
|     o pointer to refs_samRef struct to hold id/length
|     o is sorted by reference id (use functions in
|       ../genLib/strAry to access ids)
|   - samSTPtr:
|     o pointer to samEntry struct to hold line after
|       last reference (used for reading sam file)
|   - samFILE:
|     o FILE pointer to sam file header
|   - outFILE:
|     o FILE pointer to print all headers to (0 no print)
|   - headStrPtr:
|     o pointer to c-string to hold non-reference headers
|     o use null for if not saving headers
|   - lenHeadULPtr:
|     o unsigned long with headStrPtr length
| Output:
|   - Modifies:
|     o refSTPtr to have reference length and ids
|     o samSTPtr to have first read
|     o samFILE to point to first read after header
|     o outFILE to have header (or not use if 0)
|     o headStrPtr to have non-reference headers
|     o lenHeadULPtr have new headStrPtr size (if resized)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_samEntry for memory errors
|     o def_fileErr_samEntry for file errors
|     o 64 if no header reference ids in the sam file
\-------------------------------------------------------*/
signed char
getRefLen_samRef(
   struct refs_samRef *refSTPtr,/*holds ref lengths*/
   struct samEntry *samSTPtr,    /*for reading sam*/
   void *samFILE,                /*sam file with lengths*/
   void *outFILE,                /*print headers to*/
   signed char **headStrPtr,     /*holds non-ref header*/
   unsigned long *lenHeadULPtr   /*length of headStrPtr*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun07 TOC:
   '   - gets reference ids and length from sam header
   '   o fun07 sec01:
   '     - variable declarations
   '   o fun07 sec02:
   '     - blank and get first line of sam file
   '   o fun07 sec03:
   '     - get reference lengths
   '   o fun07 sec04:
   '     - return results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   signed char *tmpStr = 0;
   signed char *cpStr = 0;

   unsigned long headBytesUL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec02:
   ^   - blank and get first line of sam file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! refSTPtr->idAryStr)
   { /*If: user has not setup the strucutre*/
      if( setup_samRef(refSTPtr) )
         goto memErr_fun07_sec04;
   } /*If: user has not setup the strucutre*/

   else
      blank_refs_samRef(refSTPtr);

   if(! headStrPtr) ;
      /*not saving headers*/

   else if(! *headStrPtr)
   { /*If: need memory*/
      *headStrPtr = malloc(4096 * sizeof(signed char));

      if(! *headStrPtr)
         goto memErr_fun07_sec04;

      *lenHeadULPtr = 4069;
   } /*If: need memory*/

   /*get first line of sam file*/
   errSC = get_samEntry(samSTPtr, samFILE);

   if(errSC)
   { /*If: had error*/
      if(errSC == def_memErr_samEntry)
         goto memErr_fun07_sec04;
      else
         goto fileErr_fun07_sec04;
   } /*If: had error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec03:
   ^   - get reference lengths
   ^   o fun07 sec03 sub01:
   ^     - check if have sequence header
   ^   o fun07 sec03 sub02:
   ^     - check if have enough memory for header
   ^   o fun07 sec03 sub03:
   ^     - copy reference id and move to length
   ^   o fun07 sec03 sub04:
   ^     - get length and incurment ref number
   ^   o fun07 sec03 sub05:
   ^     - copy no reference header
   ^   o fun07 sec03 sub06:
   ^     - get next entry
   ^   o fun07 sec03 sub07:
   ^     - final error check and sort ids
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun07 Sec03 Sub01:
   *   - check if have sequence header
   \*****************************************************/

   while(! errSC)
   { /*Loop: get reference entries*/
       if(samSTPtr->extraStr[0] != '@') 
          break; /*off of header*/

       if(
             samSTPtr->extraStr[1] == 'S'
          && samSTPtr->extraStr[2] == 'Q'
          && samSTPtr->extraStr[3] == '\t'
       ){ /*If: has length information*/
          /*get past "@SQ\t"*/
          tmpStr = samSTPtr->extraStr + 4; 

          while(*tmpStr++ != ':')
          { /*Loop: past next part of SQ line*/
             if(*tmpStr < 31)
                goto nextEntry_fun07_sec03_sub06;
          } /*Loop: past next part of SQ line*/

          /**********************************************\
          * Fun07 Sec03 Sub02:
          *   - check if have enough memory for header
          \**********************************************/

          if(refSTPtr->numRefUI >= refSTPtr->arySizeUI)
          { /*If: need more memory*/
             errSC =
                realloc_refs_samRef(
                   refSTPtr,
                   refSTPtr->numRefUI + 16
                );

                if(errSC)
                   goto memErr_fun07_sec04;
          } /*If: need more memory*/

          /**********************************************\
          * Fun07 Sec03 Sub03:
          *   - copy reference id and move to length
          \**********************************************/

          
          cpStr = tmpStr;

          while(*tmpStr != '\t' && *tmpStr != '\0')
             ++tmpStr;

          errSC = *tmpStr;
          *tmpStr = '\0';

          add_strAry(
             cpStr,
             refSTPtr->idAryStr,
             refSTPtr->numRefUI
          ); /*sorting at end for better speed*/

          *tmpStr = errSC;
          ++tmpStr;
         
          if(*tmpStr < 31 )
             goto fileErr_fun07_sec04; /*end of line*/

          /*Move past the LN: flag*/
          while(*tmpStr++ !=':')
          { /*Loop: get past LN: flag*/
             if(*tmpStr < 31)
                goto fileErr_fun07_sec04; /*non-numeric*/
          } /*Loop: get past LN: flag*/

          /**********************************************\
          * Fun07 Sec03 Sub04:
          *   - get length and incurment ref number
          \**********************************************/

          /*Get the reference length*/
          tmpStr +=
             strToUI_base10str(
                tmpStr,
                &refSTPtr->lenAryUI[refSTPtr->numRefUI]
             );

          if(*tmpStr > 31)
             goto fileErr_fun07_sec04; /*non-numeric*/

          ++refSTPtr->numRefUI; /*will end at index 1*/
       } /*If: has length information*/

       /*************************************************\
       * Fun07 Sec03 Sub05:
       *   - copy no reference header
       \*************************************************/

       else if(headStrPtr)
       { /*Else: is non-reference header*/

          if(
                headBytesUL + samSTPtr->extraLenUI + 1
             >= *lenHeadULPtr
          ){ /*If: need to resize*/
             *lenHeadULPtr <<= 1;
             tmpStr =
                realloc(
                   *headStrPtr,
                     ((*lenHeadULPtr << 1) + 1)
                   * sizeof(signed char)
                );

             if(! tmpStr)
                goto memErr_fun07_sec04;

            *headStrPtr = tmpStr;
            *tmpStr = '\0';
            *lenHeadULPtr <<= 1;
          } /*If: need to resize*/

          cpStr = *headStrPtr + headBytesUL;

          cpLen_ulCp(
             cpStr,
             samSTPtr->extraStr,
             samSTPtr->extraLenUI
          );

          cpStr += samSTPtr->extraLenUI;
          *cpStr++ = str_endLine[0];

          if(str_endLine[1])
             *cpStr++ = str_endLine[1];

          *cpStr = '\0';

          headBytesUL += samSTPtr->extraLenUI + 1;
       } /*Else: is non-reference header*/

       /*************************************************\
       * Fun07 Sec03 Sub06:
       *   - get next entry
       \*************************************************/

       nextEntry_fun07_sec03_sub06:;

       if(outFILE)
          p_samEntry(samSTPtr, 0, outFILE);

       errSC = get_samEntry(samSTPtr, samFILE);
   } /*Loop: get reference entries*/

   /*****************************************************\
   * Fun07 Sec03 Sub07:
   *   - final error check and sort ids
   \*****************************************************/

   if(errSC == def_memErr_samEntry)
      goto memErr_fun07_sec04;
   else if(! refSTPtr->numRefUI)
      goto noRef_fun07_sec04;

   sortSync_strAry(
      refSTPtr->idAryStr,
      refSTPtr->lenAryUI,
      refSTPtr->numRefUI
   ); /*sort by reference id name (quicker find later)*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun07 Sec04:
   ^   - return results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun07_sec04;

   memErr_fun07_sec04:;
      errSC = def_memErr_samEntry;
      goto ret_fun07_sec04;

   fileErr_fun07_sec04:;
      errSC = def_fileErr_samEntry;
      goto ret_fun07_sec04;

   noRef_fun07_sec04:;
      errSC = 64;
      goto ret_fun07_sec04;

   ret_fun07_sec04:;
      return errSC;
} /*getRefLen_samRef*/

/*-------------------------------------------------------\
| Fun08: findRef_refs_samRef
|   - finds a reference id in a refs_samRef struct
| Input:
|   - idStr:
|     o c-string with reference id to find
|   - refSTPtr:
|     o pointer to refs_samRef struct with references
| Output:
|   - Returns:
|     o index of reference id if found
|     o < 0 if reference id not in list
\-------------------------------------------------------*/
signed long
findRef_refs_samRef(
   signed char *idStr,            /*id to find*/
   struct refs_samRef *refSTPtr /*holds ref lengths*/
){
   return
      find_strAry(
         refSTPtr->idAryStr,
         idStr,
         refSTPtr->numRefUI
      );
} /*findRef_refs_samRef*/

/*-------------------------------------------------------\
| Fun09: addRef_samRef
|   - adds reference information to array in refStack
| Input:
|   - idStr:
|     o c-string with id to add
|   - lenUI:
|     o length of reference sequence
|   - refsPtr:
|     o pointer to refs_samRef struct to add ref to
|   - errSCPtr:
|     o pointer to signed char to hold errors
| Output:
|   - Modifies:
|     o idAryStr in refsPtr to have idStr
|     o lenAryUI in refsPtr to have lenUI
|     o numRefUI in refsPtr to be resized if realloc used
|     o arrySizeUI in refsPtr to be incurmented by 1
|     o errSCPtr to be
|       * 0 for no error
|       * def_expand_samEntry if needed to realloc
|       * def_memErr_samEntry for memory error
|   - Returns
|     o index of reference
|     o -1 for errors
\-------------------------------------------------------*/
signed long
addRef_samRef(
   signed char *idStr,
   unsigned int lenUI,
   struct refs_samRef *refsPtr,
   signed char *errSCPtr
){
   unsigned long retUL = 0;

   if(! refSTPtr->idAryStr)
   { /*If: user has not setup the strucutre*/
      if( setup_samRef(refSTPtr) )
         goto memErr_fun09;
   } /*If: user has not setup the strucutre*/

   *errSCPtr = 0;


   if(refsPtr->numRefUI >= refsPtr->arySizeUI)
   { /*If: need more memory*/
      *errSCPtr =
         realloc_refs_samRef(
            refsPtr,
            refsPtr->arySizeUI + 16
         );

      if(*errSCPtr)
         goto memErr_fun09;

      *errSCPtr = def_expand_samEntry;
   } /*If: need more memory*/


   retUL =
      addSort_strAry(
         idStr,
         refsPtr->idAryStr,
         refsPtr->numRefUI
      );

   refsPtr->lenAryUI[retUL] = lenUI;
   ++refsPtr->numRefUI;

   return *errSCPtr;

   memErr_fun09:;
      return def_memErr_samEntry;
} /*addRef_samRef*/

/*-------------------------------------------------------\
| Fun10: buildRefMergeIndex_samRef
|   - looks for refseq ids that might be contigs from the
|     sam assembly and builds an index for each unique id
| Input:
|   - refsSTPtr:
|     o refs_samRef struct with reference ids to scan for
|       refseq ids from the same assembly
|     o needs to be sorted by name (default state)
| Output:
|   - Returns:
|     o signed int array to with the index assigned to
|       each reference or -1 if no other references group
|       with this reference
|     o 0 for memory errors
\-------------------------------------------------------*/
signed int *
buildRefMergeIndex_samRef(
   struct refs_samRef *refsSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun10 TOC:
   '   - looks for refseq ids that might be contigs from
   '     the sam assembly and builds an index for each
   '     unique id
   '   o fun10 sec01:
   '     - variable declarations
   '   o fun10 sec02:
   '     - initialize memory and setup pointers
   '   o fun10 sec03:
   '     - copy first id
   '   o fun10 sec04:
   '     - scan for refseq ids with multiple contigs
   '   o fun10 sec05:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *idStr = 0;
   signed char *oldIdStr = 0;
   signed char *idMemHeapStr = 0;
   unsigned int uiId = 0;
   signed int *retHeapArySI = 0;
   signed int contigLenSI = 0;
   signed int refIdSI = 0;

   signed int lenSI = 0;

   unsigned long *idULPtr = 0;
   unsigned long *oldIdULPtr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec02:
   ^   - initialize memory and setup pointers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retHeapArySI =
      malloc(refsSTPtr->numRefUI * sizeof(signed int));
   if(! retHeapArySI)
      goto memErr_fun10_sec05;

   idMemHeapStr = malloc(41 * sizeof(signed char));
      /*41 = 32 + 10; worst case is 15 offset for first
      `  id, so 15 + 16 = 31, which  means to get 10
      `  characters I need 31 + 10 characters to have
      `  two 16 byte aligned 10 char strings
      */
   if(! idMemHeapStr)
      goto memErr_fun10_sec05;

   /*make sure the id c-strings are 16 byte aligned*/
   uiId = ( ((unsigned long) idMemHeapStr) % 16 );
   if(uiId) uiId += (16 - uiId);

   idStr = idMemHeapStr + uiId;
   oldIdStr = idStr + 16;
   idStr[9] = 0;
   oldIdStr[9] = 0;

   /*assign id strings to longs for quick comparisons*/
   idULPtr = (unsigned long *) idStr;
   oldIdULPtr = (unsigned long *) oldIdStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec03:
   ^   - copy first id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   lenSI = endStr_ulCp(refsSTPtr->idAryStr);

   refIdSI = 0;

   /*copy first id if it is long enough to have contigs*/
   if(lenSI <= 13)
   { /*If: to short to be a refseq id*/
      oldIdStr[0] = 0;
      refIdSI = 1;
      retHeapArySI[0] = -1;
   } /*If: to short to be a refseq id*/

   else
   { /*Else: could be a refseq id*/
      cpLen_ulCp(
         oldIdStr,
         get_strAry(refsSTPtr->idAryStr, 0),
         9 /*assembly is NZ_<6_digits><contig_number>*/
      );

      if(oldIdStr[2] != '_')
      { /*If: not a refseq id*/
         oldIdStr[0] = 0;
         retHeapArySI[0] = -1;
         refIdSI = 1;
      } /*If: not a refseq id*/

      else
         retHeapArySI[0] = refIdSI;
   } /*Else: could be a refseq id*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec04:
   ^   - scan for refseq ids with multiple contigs
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*start scan for matching contigs*/
   for(uiId = 1; uiId < refsSTPtr->numRefUI; ++uiId)
   { /*Loop: scan for refseq ids*/
      lenSI =
        endStr_ulCp(get_strAry(refsSTPtr->idAryStr,uiId));

      if(lenSI <= 13)
      { /*If: not refseq or non-assembly refseq*/
         if(! contigLenSI)
            retHeapArySI[uiId - 1] = -1;
            /*previous sequence had one contig*/

         retHeapArySI[uiId] = -1;
         oldIdStr[0] = 0;
         refIdSI = uiId + 1;
      } /*If: not refseq or non-assembly refseq*/

      else
      { /*Else: have a refseq id*/
         cpLen_ulCp(
            idStr,
            get_strAry(refsSTPtr->idAryStr, uiId),
            9
         ); /*copy the id (<tag_2_letters>_<id_6_chars>)
            `  of the refseq id
            `  using ulCp becuase I do not know if the
            `    strings are 8byte aligned
            */

         if(
               *idULPtr != *oldIdULPtr
            || idStr[8] != oldIdStr[8]
         ){ /*If: different ids*/
            if(! contigLenSI)
               retHeapArySI[uiId - 1] = -1;
               /*previous sequence had one contig*/
            else
               refIdSI = uiId;

            retHeapArySI[uiId] = refIdSI;
            contigLenSI = 0;
            *oldIdULPtr = *idULPtr;
            oldIdStr[8] = idStr[8];
         } /*If: different ids*/

         else if(idStr[2] != '_')
         { /*Else If: not refseq id*/
            if(! contigLenSI)
               retHeapArySI[uiId - 1] = -1;
               /*previous sequence had one contig*/

            contigLenSI = 0;
            retHeapArySI[uiId] = -1;
            oldIdStr[0] = idStr[0];
            refIdSI = uiId + 1;
         } /*Else If: not refseq id*/

         else
         { /*Else: both sequences are from same assembly*/
            retHeapArySI[uiId] = refIdSI;
            ++contigLenSI;
         } /*Else: both sequences are from same assembly*/
      } /*Else: have a refseq id*/
   } /*Loop: scan for refseq ids*/

   if(! contigLenSI)
      retHeapArySI[uiId - 1] = -1;
      /*last sequence had one contig*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun10 Sec05:
   ^   - clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun10_sec05;

   memErr_fun10_sec05:;
      if(retHeapArySI)
          free(retHeapArySI);
      retHeapArySI = 0;
      goto ret_fun10_sec05;
   ret_fun10_sec05:;
      if(idMemHeapStr)
         free(idMemHeapStr);
      idMemHeapStr = 0;
      return retHeapArySI;
} /*buildRefMergeIndex_samRef*/

/*-------------------------------------------------------\
| Fun11: pSamHeader_samRef
|   - prints samEntry header for a reference or set of
|     references
| Input:
|   - refIndexSI:
|     o reference printing header to
|   - headStr:
|     o c-string with general header to print out
|     o these are headers that do not relate to the
|       reference
|     o this is made by getRefLen_samRef
|   - unmapBl:
|     o 1: print the unmapped reads header (no sequence
|          (reference) entry)
|     o 0: print header for a mapped read
|   - refSTPtr:
|     o samRef struct pointer with reference ids
|     o 0/null, use refIdStr in samSTPtr
|   - refIndexArySI:
|     o index array with reference index to print to
|     o no input (0/null), use refIdStr in samSTPtr
|     o index values:
|       * -1: no merging
|       * >=0: print reference to reference id at this
|         index
|     o must be in same order as index's in refSTPtr and
|       mereged references must be next to each other
|   - outFILE:
|     o FILE pionter to file to print to
| Output:
|   - Prints:
|     o header and reference lengths to outFILE
\-------------------------------------------------------*/
void
pSamHeader_samRef(
   signed int refIndexSI,     /*check if file was opened*/
   signed char *headStr,      /*header to print out*/
   signed char unmapBl,       /*1: unmapped header*/
   struct refs_samRef *refSTPtr, /*references used*/
   signed int *refIndexArySI, /*index to print read to*/
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun11 TOC:
   '   - prints samEntry header for a reference or set of
   '     references
   '   o fun11 sec01:
   '     - variable declarations
   '   o fun11 sec02:
   '     - print HD entry in the general header
   '   o fun11 sec03:
   '     - print sequence (reference) entries for header
   '   o fun11 sec04:
   '     - finish printing the general header
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *tmpHeadStr = headStr;
   signed int siStart = 0;
   signed int mainRefSI = 0;
   signed char tmpSC = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec02:
   ^   - print HD entry in the general header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(headStr)
   { /*If: have header to print*/
      /*as a rule the HD entry comes first, then the
      `  reference lengths (for minimap2)
      */
      if(
            tmpHeadStr[1] == 'H'
         && tmpHeadStr[2] == 'D'
         && tmpHeadStr[3] < 33
      ){ /*If: have an HD header entry*/
         tmpHeadStr += endLine_ulCp(tmpHeadStr);
         tmpSC = *tmpHeadStr;
         *tmpHeadStr = 0;

         fprintf(
           (FILE *) outFILE,
           "%s%s",
           headStr,
           str_endLine
         );
         *tmpHeadStr = tmpSC;

         while(*tmpHeadStr && *tmpHeadStr < 33)
            ++tmpHeadStr; /*get of \n and \r*/
         headStr = tmpHeadStr;
      }  /*If: have an HD header entry*/
   } /*If: have header to print*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec03:
   ^   - print sequence (reference) entries for header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(unmapBl)
      ; /*print unmapped reads; do not print seq header*/
   else if(refIndexArySI && refIndexArySI[refIndexSI]>=0)
   { /*If: user is grouping references*/
      mainRefSI = refIndexArySI[refIndexSI];

      for(siStart=refIndexSI - 1; siStart >= 0; --siStart)
         if(refIndexArySI[siStart] != mainRefSI)
            break; /*find first reference in output file*/
      ++siStart;

      while(siStart < (signed int) refSTPtr->numRefUI)
      { /*Loop: print reference headers*/
         fprintf(
            (FILE *) outFILE,
            "@SQ\tSN:%s\tLN:%u%s",
            get_strAry(refSTPtr->idAryStr, siStart),
            refSTPtr->lenAryUI[siStart],
            str_endLine
         );
         ++siStart;

         if(siStart >= (signed int) refSTPtr->numRefUI)
            break;
         if(refIndexArySI[siStart] != mainRefSI)
            break;
      } /*Loop: print reference headers*/
   } /*If: user is grouping references*/

   else if(! unmapBl)
   { /*Else: only printing this reference to the header*/
      fprintf(
         (FILE *) outFILE,
         "@SQ\tSN:%s\tLN:%u%s",
         get_strAry(refSTPtr->idAryStr, refIndexSI),
         refSTPtr->lenAryUI[refIndexSI],
         str_endLine
      );
   } /*Else: only printing this reference to the header*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun11 Sec04:
   ^   - finish printing the general header
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(headStr)
   { /*If: have general header to finish printing*/
      while(*headStr)
      { /*Loop: finish printing the general header*/
         tmpHeadStr += endLine_ulCp(tmpHeadStr);
         tmpSC = *tmpHeadStr;
         *tmpHeadStr = 0;

         fprintf(
           (FILE *) outFILE,
           "%s%s",
           headStr,
           str_endLine
         );

         *tmpHeadStr = tmpSC;
         while(*tmpHeadStr && *tmpHeadStr < 33)
            ++tmpHeadStr; /*get of \n and \r*/
         headStr = tmpHeadStr;
      } /*Loop: finish printing the general header*/
   } /*If: have general header to finish printing*/
} /*pSamHeader_samRef*/

/*-------------------------------------------------------\
| Fun12: pSamToRef_samRef
|   - print samEntry to its reference bin
| Input:
|   - samSTPtr:
|     o samEntry struct pointer with read to print
|   - prefixStr:
|     o c-string with prefix to print to
|   - headStr:
|     o c-string with header to print to the new file
|   - pUnmapBl:
|     o 1: print unmapped reads
|     o 0: do not print unmapped reads
|   - refSTPtr:
|     o samRef struct pointer with reference ids
|     o 0/null, use refIdStr in samSTPtr
|   - refIndexArySI:
|     o index array with reference index to print to
|     o no input (0/null), use refIdStr in samSTPtr
|     o index values:
|       * -1: no merging
|       * >=0: print reference to reference id at this
|         index
| Output:
|   - Prints:
|     o read in samSTPtr to file named
|       "prefix-referenceId.sam"
|   - Returns:
|     o 0 for success
|     o 1 if no reference id
|     o 2 if no sequence
|     o -1 for file errors
\-------------------------------------------------------*/
signed char
pSamToRef_samRef(
   struct samEntry *samSTPtr,  /*read to print*/
   signed char *prefixStr,     /*file prefix to print to*/
   signed char *headStr,       /*header for new file*/
   signed char pUnmapBl,       /*1: print unmapped reads*/
   struct refs_samRef *refSTPtr,/*has reference ids*/
   signed int *refIndexArySI   /*index to print read to*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun12 TOC
   '   - print samEntry to its reference bin
   '   o fun12 sec01:
   '     - variable declarations
   '   o fun12 sec02:
   '     - find (build) the file name to print the read to
   '   o fun12 sec03:
   '     - open file and print read
   '   o fun12 sec04:
   '     - return the result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char outStr[4096];
   signed int lenSI = 0;
   signed int indexSI = 0;

   signed char unmapBl = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec02:
   ^   - find (build) the file name to print the read to
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   lenSI = cpStr_ulCp(outStr, prefixStr);
   outStr[lenSI++] = '-';
   outStr[lenSI] = 0;

   if(! samSTPtr->seqStr[0] || samSTPtr->seqStr[0] == '*')
      goto noSeq_fun12;

   else if(pUnmapBl && samSTPtr->flagUS & 4)
   { /*Else If: unmapped read and print unmapped reads*/
      lenSI +=
          cpStr_ulCp(
             &outStr[lenSI],
             (signed char *) "unmapped"
          );
      unmapBl = 1;
   } /*Else If: unmapped read and print unmapped reads*/

   else if(! refIndexArySI || ! refSTPtr)
      lenSI +=
          cpStr_ulCp(&outStr[lenSI], samSTPtr->refIdStr);
   else
   { /*Else: find id in array*/
      indexSI =
         find_strAry(
            refSTPtr->idAryStr,
            samSTPtr->refIdStr,
            refSTPtr->numRefUI
         );
      if(indexSI < 0)
         goto noId_fun12;
      else if(refIndexArySI[indexSI] >= 0)
         indexSI = refIndexArySI[indexSI];

      lenSI +=
          cpStr_ulCp(
             &outStr[lenSI],
             get_strAry(refSTPtr->idAryStr, indexSI)
         );
   } /*Else: find id in array*/

   outStr[lenSI++] = '.';
   outStr[lenSI++] = 's';
   outStr[lenSI++] = 'a';
   outStr[lenSI++] = 'm';
   outStr[lenSI] = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec03:
   ^   - open file and print read
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   outFILE = fopen((char *) outStr, "r");

   if(! outFILE)
   { /*If: need to make the file*/
      outFILE = fopen((char *) outStr, "w");
      if(outFILE)
         pSamHeader_samRef(
            indexSI,
            headStr,
            unmapBl,
            refSTPtr,
            refIndexArySI,
            outFILE
         ); /*print header to sam file*/
   } /*If: need to make the file*/

   else
   { /*Else: file already exists*/
      fclose(outFILE);
      outFILE = fopen((char *) outStr, "a");
   } /*Else: file already exists*/

   if(! outFILE)
      goto fileErr_fun12;
   p_samEntry(samSTPtr, 0, outFILE);
   fclose(outFILE);
   outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec04:
   ^   - return the result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   noId_fun12:;
      return 1;
   noSeq_fun12:;
      return 2;
   fileErr_fun12:;
      return -1;
} /*pSamToRef_samRef*/

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
