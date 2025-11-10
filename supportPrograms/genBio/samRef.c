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
'   o fun07: getRefLen_samEntry
'     - gets reference ids & length from a sam file header
'   o fun08: findRef_refs_samRef
'     - finds a reference id in a refs_samRef struct
'   o fun09: addRef_samRef
'     - adds reference information to array in refStack
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
| Fun07: getRefLen_samEntry
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
\-------------------------------------------------------*/
signed char
getRefLen_samEntry(
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

   ret_fun07_sec04:;
   return errSC;
} /*getRefLen_samEntry*/

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
