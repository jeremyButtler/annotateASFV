/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' checkSum TOC:
'   - library for checksum computing
'   o header:
'     - defiend variables and guards
'   o .c tbl01: crc32Tbl_checkSum
'     - crc32 lookup table
'     - copied from misc0110's libattpng repository
'   o fun01: usAdler32_checkSum
'     - find adler32 checksum with unsigned shorts; 16 bit
'   o fun02: usAdler32Byte_checkSum
'     - does adler32 (signed short) checksum on one byte
'   o fun03: usAdler32Finish_checkSum
'     - does little endin swap for final value of adler32
'   o fun04: adler32_checkSum
'     - finds adler32 checksum with ints; 32 bit
'   o fun05: adler32Byte_checkSum
'     - adds a byte to an adler32 checksum
'   o .h fun06: adler32Finish_checkSum
'     - finishes the addler32 checksum
'   o .c fun07: crc32Mac_checkSum
'     - find new crc value after adding next byte
'     - copied from misc0110's libattpng repository, I
'       really do not know whats going on here
'   o fun08: crc32_checkSum
'     - does a crc32 checksum on input data
'   o fun09: crc32Byte_checkSum
'     - adds byte to a crc32 checksum
'   o .h fun10: crc32Finish_checkSum
'     - does final step in crc32 checkusm
'   o fun11: mkCrc32Tbl_checkSum
'     - makes a crc32 table (used to generate this table)
'     - this was adapted from python code from github
'       co-pilot, which it in turn took from others
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - defiend variables and guards
\-------------------------------------------------------*/

#ifndef CHECKSUM_H
#define CHECKSUM_H

#define def_adler32_checkSum 65521
#define def_adler32ToUSMax_checkSum 14
   /*max unsigned short is 65535, so + 14 would get
   `  the maximum for an unsigned short
   */

/*-------------------------------------------------------\
| Fun01: usAdler32_checkSum
|   - finds adler32 checksum with unsigned shorts; 16 bit
| Input:
|   - dataAryUC:
|     o unsigned char array with data to get adler32 with
|   - lenUL:
|     o length to find adler32 for in dataAryUC
|   - adlerAryUC:
|     o 4 byte array to get adler32 value
| Output:
|   - Modifies:
|     o adlerAryUC to have adler32 sum (little endin)
\-------------------------------------------------------*/
void
usAdler32_checkSum(
   unsigned char *dataAryUC,  /*data to sum*/
   unsigned long lenUL,       /*amount of data*/
   unsigned char *adlerAryUC  /*4 bytes, gets adler32*/
);

/*-------------------------------------------------------\
| Fun02: usAdler32Byte_checkSum
|   - does adler32 (signed short) checksum on one byte
| Input:
|   - byteUC:
|     o unsigned char to add to adler32 checksum
|   - adlerAryUC:
|     o 4 byte array to update adler32 checksum
| Output:
|   - Modifies:
|     o adlerAryUC to have updated adler32 sum
| Note:
|  - you will need to convert to little endin at end
|    unsigned short *oneUS = (unsigned short *) &aryUC[2];
|    unsigned short *twoUS = (unsigned short *) &aryUC[0];
|    *oneUS = usToLittle_endin(*oneUS);
|    *twoUS = usToLittle_endin(*twoUS);
\-------------------------------------------------------*/
void
usAdler32Byte_checkSum(
   unsigned char byteUC,      /*byte to add*/
   unsigned char *adlerAryUC  /*4 bytes, gets adler32*/
);

/*-------------------------------------------------------\
| Fun03: usAdler32Finish_checkSum
|   - does little endin swap for final value of adler32
| Input:
|   - adlerAryUC:
|     o 4 byte array with adler32 checksum to finish
| Output:
|   - Modifies:
|     o adlerAryUC to have final adler32 format
\-------------------------------------------------------*/
void
usAdler32Finish_checkSum(
   unsigned char *aryUC    /*4 bytes, value to finish*/
);

/*-------------------------------------------------------\
| Fun04: adler32_checkSum
|   - finds adler32 checksum with ints; 32 bit
| Input:
|   - dataAryUC:
|     o unsigned char array with data to get adler32 with
|   - lenUL:
|     o length to find adler32 for in dataAryUC
| Output:
|   - Returns:
|     o adler32 calculated checksum
\-------------------------------------------------------*/
signed int
adler32_checkSum(
   unsigned char *dataAryUC,  /*data to sum*/
   unsigned long lenUL        /*amount of data*/
);

/*-------------------------------------------------------\
| Fun05: adler32Byte_checkSum
|   - adds a byte to an adler32 checksum
| Input:
|   - byteUC:
|     o byte to add to adler32 checksum
|   - oneSIPtr:
|     o signed int pointer with first checksum value
|   - twoSIPtr:
|     o signed int pointer with second checksum value
| Output:
|   - Modifies:
|     o oneSIPtr to have updated checksum for first value
|     o twoSIPtr to have updated checksum for 2nd value
\-------------------------------------------------------*/
void
adler32Byte_checkSum(
   unsigned char byteUC,  /*byte to add to checksum*/
   signed int *oneSIPtr,  /*first value in adler32*/
   signed int *twoSIPtr   /*second value in adler32*/
);

/*-------------------------------------------------------\
| Fun06: adler32Finish_checkSum
|   - finishes the addler32 checksum
| Input:
|   - oneMac:
|     o first value in adler32 checksum (signed int)
|   - twoSIPtr:
|     o second value in adler32 checksum (signed int)
| Output:
|   - Returns:
|     o adler32 checksum
\-------------------------------------------------------*/
#define adler32Finish_checkSum(oneMac, twoMac) ( ((twoMac) << 16) | (oneMac))

/*-------------------------------------------------------\
| Fun08: crc32_checkSum
|   - does a crc32 checksum on input data
| Input:
|   - aryUC:
|     o unsigned char array with data to get crc32 with
|   - lenUL:
|     o length of dataAryUC
| Output:
|   - Returns:
|     o crc32 checksum
\-------------------------------------------------------*/
unsigned int
crc32_checkSum(
   unsigned char *aryUC, /*data to get crc32 checksum on*/
   unsigned long lenUL   /*bytes in aryUC*/
);

/*-------------------------------------------------------\
| Fun09: crc32Byte_checkSum
|   - adds byte to a crc32 checksum
| Input:
|   - byteUC:
|     o byte to add to crc32 checksum
|   - crc32UI:
|     o currnet crc32 checksum
| Output:
|   - Returns:
|     o updated crc32 checksum
\-------------------------------------------------------*/
unsigned int
crc32Byte_checkSum(
   unsigned char byteUC, /*data to add to crc32 checksum*/
   unsigned int crc32UI  /*current crc32 checksum*/
);

/*-------------------------------------------------------\
| Fun10: crc32Finish_checkSum
|   - does final step in crc32 checkusm
| Input:
|   - crc32Mac:
|     o crc32 checksum to finalize
| Output:
|   - Returns:
|     o finalized crc32 checksum
\-------------------------------------------------------*/
#define crc32Finish_checkSum(crc32Mac) ((unsigned int) ~(crc32Mac))

/*-------------------------------------------------------\
| Fun11: mkCrc32Tbl_checkSum
|   - makes a crc32 table (used to generate this table)
|   - this was adapted from python code from github
|     co-pilot, which it in turn took from others
| Input:
|   - crc32AryUI:
|     o unsigned int array to hold crc32 values; needs to
|       have at least 256 elements
| Output:
|   - Modifies:
|     o crc32AryUI to have the crc32 table
\-------------------------------------------------------*/
void
mkCrc32Tbl_checkSum(
   unsigned int *crc32AryUI
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
