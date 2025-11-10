/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' bigEndin SOF: Start Of File
'   - converts 
'   o header:
'     - guards
'   o fun01: ulToBig_endin
'     - converts unsigned long to big endin format
'   o fun02: uiToBig_endin
'     - converts unsigned int to big endin format
'   o fun03: usToBig_endin
'     - converts unsigned short to big endin format
'   o fun04: ulToLittle_endin
'     - converts unsigned long to little endin format
'   o fun05: uiToLittle_endin
'     - converts unsigned int to little endin format
'   o fun06: usToLittle_endin
'     - converts unsigned short to little endin format
'   o fun07: ulSwap_endin
'     - converts unsigned long to different endin format
'   o fun08: uiSwap_endin
'     - converts unsigned int to oppisite endin format
'   o fun09: usSwap_endin
'     - converts unsigned short to oppisite endin format
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef BIG_ENDIN_CONVERSION_H
#define BIG_ENDIN_CONVERSION_H

/*-------------------------------------------------------\
| Fun01: ulToBig_endin
|   - converts unsigned long to big endin format
| Input:
|   - valUL:
|     o unsigned long to convert
| Output:
|   - Returns:
|     o bigendin number as unsigned long
\-------------------------------------------------------*/
unsigned long
ulToBig_endin(
   unsigned long valUL
);

/*-------------------------------------------------------\
| Fun02: uiToBig_endin
|   - converts unsigned int to big endin format
| Input:
|   - valUI:
|     o unsigned int to convert
| Output:
|   - Returns:
|     o bigendin number as unsigned int
\-------------------------------------------------------*/
unsigned int
uiToBig_endin(
   unsigned int valUI
);

/*-------------------------------------------------------\
| Fun03: usToBig_endin
|   - converts unsigned short to big endin format
| Input:
|   - valUI:
|     o unsigned short to convert
| Output:
|   - Returns:
|     o bigendin number as unsigned short
\-------------------------------------------------------*/
unsigned short
usToBig_endin(
   unsigned short valUS
);

/*-------------------------------------------------------\
| Fun04: ulToLittle_endin
|   - converts unsigned long to little endin format
| Input:
|   - valUL:
|     o unsigned long to convert
| Output:
|   - Returns:
|     o little endin number as unsigned long
\-------------------------------------------------------*/
unsigned long
ulToLittle_endin(
   unsigned long valUL
);

/*-------------------------------------------------------\
| Fun05: uiToLittle_endin
|   - converts unsigned int to little endin format
| Input:
|   - valUI:
|     o unsigned int to convert
| Output:
|   - Returns:
|     o little endin number as unsigned int
\-------------------------------------------------------*/
unsigned int
uiToLittle_endin(
   unsigned int valUI
);

/*-------------------------------------------------------\
| Fun06: usToLittle_endin
|   - converts unsigned short to little endin format
| Input:
|   - valUI:
|     o unsigned short to convert
| Output:
|   - Returns:
|     o little endin number as unsigned short
\-------------------------------------------------------*/
unsigned short
usToLittle_endin(
   unsigned short valUS
);

/*-------------------------------------------------------\
| Fun07: ulSwap_endin
|   - converts unsigned long to different endin format
| Input:
|   - valUL:
|     o unsigned long to convert
| Output:
|   - Returns:
|     o oppisite endin format as valUL
\-------------------------------------------------------*/
unsigned long
ulSwap_endin(
   unsigned long valUL
);

/*-------------------------------------------------------\
| Fun08: uiSwap_endin
|   - converts unsigned int to oppisite endin format
| Input:
|   - valUI:
|     o unsigned int to convert
| Output:
|   - Returns:
|     o oppisiste endin number as unsigned int
\-------------------------------------------------------*/
unsigned int
uiSwap_endin(
   unsigned int valUI
);

/*-------------------------------------------------------\
| Fun09: usSwap_endin
|   - converts unsigned short to oppisite endin format
| Input:
|   - valUI:
|     o unsigned short to convert
| Output:
|   - Returns:
|     o oppisite endin number as unsigned short
\-------------------------------------------------------*/
unsigned short
usSwap_endin(
   unsigned short valUS
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
