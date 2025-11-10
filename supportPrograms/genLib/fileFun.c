/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' fileFun SOF: Start Of File
'   - holds functions for dealing with non-gz files
'   o header:
'     - included libraries
'   o fun01: lineCnt_fileFun
'     - count number of lines in a file
'   o fun02: getLine_fileFun
'     - get line or x number of bytes from file
'     - this may be slower than fgets, but it also avoids
'       line break issues between OS's
'   o fun03: getFullLine_fileFun
'     - gets a full length line from a file
'   o license:
'     - licensing for this code (public domain / mit)
'   o options:
'     - -DSLOW_BREAK sets fileFun to not use fgets, but
'       instead the fgetc, which is slower, but also
'       more reliable system
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
#include "../genLib/ulCp.h"

/*.h files only*/
#include "endLine.h"

/*-------------------------------------------------------\
| Fun01: lineCnt_fileFun
|   - count number of lines in a file
| Input:
|   - inFILE:
|     o FILE pointer to get number of lines in
|   - lenSLPtr:
|     o signed long pointer to hold length of longest line
| Output:
|   - Modifies:
|     o inFILE to be at start of file (SEEK_SET)
|     o lenSLPtr to have length of longest line
|       * does not include line break (`\n`, `\r\n`, `\r`,
|         or `\n\r`)
|   - Returns:
|     o number of lines (index 1) in file
|       * EOF is counted if file does not end in a line
|         break
\-------------------------------------------------------*/
signed long
lineCnt_fileFun(
   void *inFILE,         /*get number lines in inFILE*/
   signed long *lenSLPtr /*length of longest line*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - counts number of lines in file
   '   o fun01 sec01:
   '     - variable declarations
   '   o fun01 sec02:
   '     - count number of lines
   '   o fun01 sec03:
   '     - finalize break counts and max length & return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define size_fun01 1024
   signed char buffStr[size_fun01 + 8];
   unsigned int tmpUI = 0;
   signed long posSL = 0;
   signed long lenSL = 0;
   signed long lineCntSL = 0;
   signed char noEofBreakBl = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - count number of lines
   ^   o fun01 sec02 sub01:
   ^     - read in first part of file
   ^   o fun01 sec02 sub02:
   ^     - find line break in buffer + start count loop
   ^   o fun01 sec02 sub03:
   ^     - max length & count breaks starting with \n
   ^   o fun01 sec02 sub04:
   ^     - count breaks starting with carriage returns
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *   - read in first part of file
   \*****************************************************/

   *lenSLPtr = 0;

   noEofBreakBl = 0;
   tmpUI =
      fread(
         (char *) buffStr,
         sizeof(signed char),
         size_fun01,
         (FILE *) inFILE
      );

   buffStr[tmpUI] = 0;
   if(! tmpUI)
      goto endFile_fun01_sec03;
   tmpUI = 0;

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *   - find next line break in buffer + start count loop
   \*****************************************************/

   while(1 == 1)
   { /*Loop: find number of lines*/
      tmpUI = endLine_ulCp(&buffStr[posSL]);
         /*always ends at \n, \r, or \0 only*/
      posSL += tmpUI;
      lenSL += tmpUI;

      if(buffStr[posSL] == '\n')
         ;
      else if(buffStr[posSL] == '\r')
         ;
      else
      { /*Else: not at line break, so end of string*/
         noEofBreakBl = 0;
         tmpUI =
            fread(
               (char *) buffStr,
               sizeof(signed char),
               size_fun01,
               (FILE *) inFILE
            );

         buffStr[tmpUI] = 0;
         if(! tmpUI)
            goto endFile_fun01_sec03;
         posSL = 0;
         continue;
      } /*Else: not at line break, so end of string*/

      /**************************************************\
      * Fun01 Sec02 Sub03:
      *     - max length & count breaks starting with \n
      \**************************************************/

      if(lenSL > *lenSLPtr)
         *lenSLPtr = lenSL;
      noEofBreakBl = 0;
      lenSL = 0;

      while(buffStr[posSL] == '\n')
      { /*Loop: count number of new line endings*/
         ++lineCntSL;
         ++posSL;

         if(buffStr[posSL] == '\r')
            ++posSL;

         if(buffStr[posSL] == '\0')
         { /*Else If: at end of buffer*/
            tmpUI =
              fread(
                 (char *) buffStr,
                 sizeof(signed char),
                 3,
                 (FILE *) inFILE
              );

            buffStr[tmpUI] = 0;
            if(! tmpUI)
               goto endFile_fun01_sec03;

            posSL = 0;

            if(buffStr[posSL] == '\r')
            { /*If: have \n\r line break*/
               ++posSL;
               continue;
            } /*If: have \n\r line break*/
         } /*Else If: at end of buffer*/
      } /*Loop: count number of new line endings*/

      /**************************************************\
      * Fun01 Sec02 Sub04:
      *   - count number breaks with carriage returns
      \**************************************************/

      while(buffStr[posSL] == '\r')
      { /*Loop: count number of new line endings*/
         ++lineCntSL;
         ++posSL;

         if(buffStr[posSL] == '\n')
            ++posSL;

         if(buffStr[posSL] == '\0')
         { /*Else If: at end of buffer*/
            tmpUI =
              fread(
                 (char *) buffStr,
                 sizeof(signed char),
                 3,
                 (FILE *) inFILE
              );

            buffStr[tmpUI] = 0;
            if(! tmpUI)
               goto endFile_fun01_sec03;

            posSL = 0;

            if(buffStr[posSL] == '\n')
            { /*If: have \r\n line break*/
               ++posSL;
               continue;
            } /*If: have \r\n line break*/
         } /*Else If: at end of buffer*/
      } /*Loop: count number of new line endings*/
   } /*Loop: find number of lines*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - finalize break counts and max length & return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   endFile_fun01_sec03:;
      lineCntSL += noEofBreakBl;
         /*account for cases were no new line at EOF*/

      if(lenSL > *lenSLPtr)
         *lenSLPtr = lenSL;

      fseek((FILE *) inFILE, 0, SEEK_SET);
      return lineCntSL;
} /*lineCnt_fileFun*/

/*-------------------------------------------------------\
| Fun02: getLine_fileFun
|   - get line or x number of bytes from file
|   - the non-quick_break method is slower than fgets, but
|     will never have an issue with a line break
| Input:
|   - inFILE:
|     o FILE pointer to get number of lines in
|   - buffStr:
|     o buffer to hold line or bytesSL of bytes
|   - bytesSL:
|     o number bytes to read (must be > 3)
|   - lenSLPtr:
|     o signed long pointer to hold number of real bytes
|       (including one to two bytes in line break) in file
| Output:
|   - Modifies:
|     o inFILE to be after next line or after bytesSL
|     o buffStr to have line or number of requested bytes
|     o lenSLPtr have the number bytes read from file
|       * this can be different then then return value,
|         which is the number bytes in buffStr includes
|   - Returns:
|     o number of bytes in buffer
|       * this will include the line break for your OS, so
|         with the -DWINDOWS flag it will always end with
|         '\r\n' even when reading a unix file (uses '\n')
|         i. or with -DUNIX, it will replace any '\r',
|            '\r\n', or '\n\r' with a '\n'.
|     o 0 for EOF
|     o -1 if bytesSL is less then 3
\-------------------------------------------------------*/
signed long
getLine_fileFun(
   void *inFILE,         /*file to get bytes from*/
   signed char *buffStr, /*buffer to get bytes*/
   signed long bytesSL,  /*number bytes to read in*/
   signed long *lenSLPtr /*actual number bytes from file*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - get line or x number of bytes from file
   '   o fun02 sec01:
   '     - variable declartions and quick checks
   '   o fun02 sec02
   '     - quick_break, fgets system
   '   o fun02 sec03:
   '     - uses fgets, which is much quicker, but also
   '       depends on fgets not missing line breaks
   '   o fun02 sec04:
   '     - return and set null
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01
   ^   - variable declartions and quick checks
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed long retSL = 0;
   *lenSLPtr = 0;

   if(bytesSL < 3)
      return -1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02
   ^   - slower getc system (will never fail)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #ifdef SLOW_BREAK
      while(retSL < bytesSL - 2)
      { /*Loop: get line*/
         buffStr[retSL] = fgetc((FILE *) inFILE);

         if(buffStr[retSL] == EOF)
            break; /*EOF*/

         else if(buffStr[retSL] == '\n')
         { /*Else If: have new line*/
            ++*lenSLPtr;
            ++retSL;
            buffStr[retSL] = fgetc((FILE *) inFILE);

            if(buffStr[retSL] == EOF)
            { /*If: hit EOF*/
               buffStr[retSL - 1] = str_endLine[0];
               if(str_endLine[1])
                  buffStr[retSL++] = str_endLine[1];
               break; /*EOF*/
            } /*If: hit EOF*/

            else if(buffStr[retSL] != '\r')
               ungetc(buffStr[retSL], (FILE *) inFILE);
            else
               ++*lenSLPtr;

            buffStr[retSL - 1] = str_endLine[0];
            if(str_endLine[1])
               buffStr[retSL++] = str_endLine[1];
            break;
         } /*Else If: have new line*/

         else if(buffStr[retSL] == '\r')
         { /*Else If: have carraige return*/
            ++*lenSLPtr;
            ++retSL;
            buffStr[retSL] = fgetc((FILE *) inFILE);

            if(buffStr[retSL] == EOF)
            { /*If: hit EOF*/
               buffStr[retSL - 1] = str_endLine[0];
               if(str_endLine[1])
                  buffStr[retSL++] = str_endLine[1];
               break; /*EOF*/
            } /*If: hit EOF*/
            
            else if(buffStr[retSL] != '\n')
               ungetc(buffStr[retSL], (FILE *) inFILE);
            else
               ++*lenSLPtr;

            buffStr[retSL - 1] = str_endLine[0];
            if(str_endLine[1])
               buffStr[retSL++] = str_endLine[1];
            break;
         } /*Else If: have carraige return*/

         else
         { /*Else: not line break*/
            ++retSL;
            ++*lenSLPtr;
         } /*Else: not line break*/
      } /*Loop: get line*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - uses fgets, which is much quicker, but also
   ^     depends on fgets not missing line breaks
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #else
      --bytesSL;

      if(fgets((char *) buffStr,bytesSL,(FILE *) inFILE))
      { /*If: have input from file*/
         retSL = endStr_ulCp(buffStr);
         *lenSLPtr= retSL;

         if(buffStr[retSL - 1] == '\n')
         { /*If: have line break*/
            --retSL;

            if(buffStr[retSL] == '\r')
               --retSL;
            else
            { /*Else: check if next character is \r*/
               buffStr[retSL] = fgetc((FILE *) inFILE);

               if(buffStr[retSL] == EOF)
                  ;
               else if(buffStr[retSL] == '\r')
                  --retSL;
               else
                  ungetc(buffStr[retSL], (FILE *) inFILE);
            } /*Else: check if next character is \r*/

            buffStr[retSL++] = str_endLine[0];
            if(str_endLine[1])
               buffStr[retSL++] = str_endLine[1];
            buffStr[retSL] = '\0';
         } /*If: have line break*/

         else if(buffStr[retSL - 1] == '\r')
         { /*If: have line break*/
            --retSL;

            if(buffStr[retSL] == '\n')
               --retSL;
            else
            { /*Else: check if next character is \n*/
               buffStr[retSL] = fgetc((FILE *) inFILE);

               if(buffStr[retSL] == EOF)
                  ;
               else if(buffStr[retSL] == '\n')
                  --retSL;
               else
                  ungetc(buffStr[retSL], (FILE *) inFILE);
            } /*Else: check if next character is \n*/

            buffStr[retSL++] = str_endLine[0];
            if(str_endLine[1])
               buffStr[retSL++] = str_endLine[1];
         } /*If: have line break*/
      } /*If: have input from file*/
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec04:
   ^   - return and set null
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   buffStr[retSL] = 0;
   return retSL;
} /*getLine_fileFun*/

/*-------------------------------------------------------\
| Fun03: getFullLine_fileFun
|   - gets a full length line from a file
| Input:
|   - inFILE:
|     o FILE pointer to grab line from
|   - buffStrPtr:
|     o c-string pointer to add buffer to (expanded if
|       needed)
|   - buffSizeSLPtr:
|     o signed long pointer with current size of buffer
|   - lenSLPtr:
|     o signed long pointer to hold number bytes read in
|       from file (not total in buffStrPtr)
|   - offsetSL:
|     o position to start at in buffStrPtr
| Output:
|   - Modifies:
|     o inFILE to be on character after line
|     o buffStrPtr to have new characters and if needed,
|       resized
|     o buffSizeSLPtr to have new size of buffStrPtr if
|       buffStrPtr is resized
|     o lenSLPtr to have number bytes read in from file
|   - Returns:
|     o number bytes in buffStrPtr
|     o 0 for EOF
|     o -1 for memory errors
\-------------------------------------------------------*/
signed long
getFullLine_fileFun(
  void *inFILE,               /*file to get output from*/
  signed char **buffStrPtr,   /*buffer to hold output*/
  signed long *buffSizeSLPtr, /*size of buffer*/
  signed long *lenSLPtr,      /*bytes read from file*/
  signed long offsetSL        /*offset in buffer*/
){
   signed char *swapStr = 0;
   signed long retSL = 0;
   signed long tmpSL = 0;
   signed long bytesSL = 0;

   *lenSLPtr = 0;

   if(! *buffStrPtr)
   { /*If: need memory*/
      *buffSizeSLPtr = 516;

      *buffStrPtr =
         calloc(*buffSizeSLPtr + 8, sizeof(signed char));
      if(! *buffStrPtr)
         goto memErr_fun03_sec0x;
   } /*If: need memory*/

	do{ /*Loop: get line*/
      bytesSL = *buffSizeSLPtr - offsetSL;

      if(bytesSL < 5)
      { /*If: need more memory*/
         tmpSL = *buffSizeSLPtr;
         tmpSL  += (*buffSizeSLPtr >> 1);
         swapStr =
            realloc(
               *buffStrPtr,
               (tmpSL + 8) * sizeof(signed char)
            );
         if(! swapStr)
            goto memErr_fun03_sec0x;
         *buffStrPtr = swapStr;
         *buffSizeSLPtr = tmpSL;
         bytesSL = *buffSizeSLPtr - offsetSL;
      } /*If: need more memory*/

      bytesSL =
         getLine_fileFun(
            inFILE,
            &(*buffStrPtr)[offsetSL],
            bytesSL,
            &tmpSL
         );
      *lenSLPtr += tmpSL;
      retSL += bytesSL;
      offsetSL += bytesSL;
   } while(
          bytesSL
       && (*buffStrPtr)[offsetSL - 1] != '\n'
       && (*buffStrPtr)[offsetSL - 1] != '\r'
     ); /*Loop: get line*/

   goto ret_fun03_sec0x;

   memErr_fun03_sec0x:;
      retSL = -1;
      goto ret_fun03_sec0x;
   ret_fun03_sec0x:;
       return retSL;
} /*getFullLine_fileFun*/

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
