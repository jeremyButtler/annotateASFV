/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' fileFun SOF: Start Of File
'   - holds functions for dealing with non-gz files
'   o header:
'     - guards
'   o fun01: lineCnt_fileFun
'     - count number of lines in a file
'   o fun02: getLine_fileFun
'     - get line or x number of bytes from file
'     - the non-quick_break method is slower than fgets,
'       but will never have an issue with a line break
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
|   - guards
\-------------------------------------------------------*/

#ifndef NON_STANDARD_FILE_FUNCTIONS
#define NON_STANDARD_FILE_FUNCTIONS

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
);

/*-------------------------------------------------------\
| Fun02: getLine_fileFun
|   - get line or x number of bytes from file
|   - the non-quick_break method is slower than fgets,
|     but will never have an issue with a line break
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
);

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
);

#endif

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
