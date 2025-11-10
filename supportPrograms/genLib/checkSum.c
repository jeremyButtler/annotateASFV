/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' checkSum TOC:
'   - library for checksum computing
'   o header:
'     - included libraries
'   o .c tbl01: crc32Tbl_checkSum
'     - crc32 lookup table
'     - copied from misc0110's libattpng repository
'   o fun01: usAdler32_checkSum
'     - find adler32 checksum with unsigned shorts; 16 bit
'     - on O3 this is about as fast as alder32_checksum,
'       but on O0 is slower
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
'   o fun07: crc32Mac_checkSum
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
' Credit:
'     - a good chunk of the crc32 from misc0110's
'       libattpng repository
'     - also https://github.com/davecom/nflate helped with
'       the crc32 sums
'     - github copilot hepled with crc32 sums
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "checkSum.h"
#include "endin.h"

/*-------------------------------------------------------\
| Tbl01: crc32Mac_checkSum
|   - crc32 lookup table
|   - copied from misc0110's libattpng repository
\-------------------------------------------------------*/
unsigned int
   crc32Tbl_checkSum[256] = {
      0x0000000, 0x77073096, 0xee0e612c, 0x990951ba,
      0x76dc419, 0x706af48f, 0xe963a535, 0x9e6495a3,
      0xedb8832, 0x79dcb8a4, 0xe0d5e91e, 0x97d2d988,
      0x09b64c2b, 0x7eb17cbd, 0xe7b82d07, 0x90bf1d91,
      0x1db71064, 0x6ab020f2, 0xf3b97148, 0x84be41de,
      0x1adad47d, 0x6ddde4eb, 0xf4d4b551, 0x83d385c7,
      0x136c9856, 0x646ba8c0, 0xfd62f97a, 0x8a65c9ec,
      0x14015c4f, 0x63066cd9, 0xfa0f3d63, 0x8d080df5,
      0x3b6e20c8, 0x4c69105e, 0xd56041e4, 0xa2677172,
      0x3c03e4d1, 0x4b04d447, 0xd20d85fd, 0xa50ab56b,
      0x35b5a8fa, 0x42b2986c, 0xdbbbc9d6, 0xacbcf940,
      0x32d86ce3, 0x45df5c75, 0xdcd60dcf, 0xabd13d59,
      0x26d930ac, 0x51de003a, 0xc8d75180, 0xbfd06116,
      0x21b4f4b5, 0x56b3c423, 0xcfba9599, 0xb8bda50f,
      0x2802b89e, 0x5f058808, 0xc60cd9b2, 0xb10be924,
      0x2f6f7c87, 0x58684c11, 0xc1611dab, 0xb6662d3d,
      0x76dc4190, 0x1db7106, 0x98d220bc, 0xefd5102a,
      0x71b18589, 0x6b6b51f, 0x9fbfe4a5, 0xe8b8d433,
      0x7807c9a2, 0xf00f934, 0x9609a88e, 0xe10e9818,
      0x7f6a0dbb, 0x86d3d2d, 0x91646c97, 0xe6635c01,
      0x6b6b51f4, 0x1c6c6162, 0x856530d8, 0xf262004e,
      0x6c0695ed, 0x1b01a57b, 0x8208f4c1, 0xf50fc457,
      0x65b0d9c6, 0x12b7e950, 0x8bbeb8ea, 0xfcb9887c,
      0x62dd1ddf, 0x15da2d49, 0x8cd37cf3, 0xfbd44c65,
      0x4db26158, 0x3ab551ce, 0xa3bc0074, 0xd4bb30e2,
      0x4adfa541, 0x3dd895d7, 0xa4d1c46d, 0xd3d6f4fb,
      0x4369e96a, 0x346ed9fc, 0xad678846, 0xda60b8d0,
      0x44042d73, 0x33031de5, 0xaa0a4c5f, 0xdd0d7cc9,
      0x5005713c, 0x270241aa, 0xbe0b1010, 0xc90c2086,
      0x5768b525, 0x206f85b3, 0xb966d409, 0xce61e49f,
      0x5edef90e, 0x29d9c998, 0xb0d09822, 0xc7d7a8b4,
      0x59b33d17, 0x2eb40d81, 0xb7bd5c3b, 0xc0ba6cad,
      0xedb88320, 0x9abfb3b6, 0x3b6e20c, 0x74b1d29a,
      0xead54739, 0x9dd277af, 0x4db2615, 0x73dc1683,
      0xe3630b12, 0x94643b84, 0xd6d6a3e, 0x7a6a5aa8,
      0xe40ecf0b, 0x9309ff9d, 0xa00ae27, 0x7d079eb1,
      0xf00f9344, 0x8708a3d2, 0x1e01f268, 0x6906c2fe,
      0xf762575d, 0x806567cb, 0x196c3671, 0x6e6b06e7,
      0xfed41b76, 0x89d32be0, 0x10da7a5a, 0x67dd4acc,
      0xf9b9df6f, 0x8ebeeff9, 0x17b7be43, 0x60b08ed5,
      0xd6d6a3e8, 0xa1d1937e, 0x38d8c2c4, 0x4fdff252,
      0xd1bb67f1, 0xa6bc5767, 0x3fb506dd, 0x48b2364b,
      0xd80d2bda, 0xaf0a1b4c, 0x36034af6, 0x41047a60,
      0xdf60efc3, 0xa867df55, 0x316e8eef, 0x4669be79,
      0xcb61b38c, 0xbc66831a, 0x256fd2a0, 0x5268e236,
      0xcc0c7795, 0xbb0b4703, 0x220216b9, 0x5505262f,
      0xc5ba3bbe, 0xb2bd0b28, 0x2bb45a92, 0x5cb36a04,
      0xc2d7ffa7, 0xb5d0cf31, 0x2cd99e8b, 0x5bdeae1d,
      0x9b64c2b0, 0xec63f226, 0x756aa39c, 0x26d930a,
      0x9c0906a9, 0xeb0e363f, 0x72076785, 0x5005713,
      0x95bf4a82, 0xe2b87a14, 0x7bb12bae, 0xcb61b38,
      0x92d28e9b, 0xe5d5be0d, 0x7cdcefb7, 0xbdbdf21,
      0x86d3d2d4, 0xf1d4e242, 0x68ddb3f8, 0x1fda836e,
      0x81be16cd, 0xf6b9265b, 0x6fb077e1, 0x18b74777,
      0x88085ae6, 0xff0f6a70, 0x66063bca, 0x11010b5c,
      0x8f659eff, 0xf862ae69, 0x616bffd3, 0x166ccf45,
      0xa00ae278, 0xd70dd2ee, 0x4e048354, 0x3903b3c2,
      0xa7672661, 0xd06016f7, 0x4969474d, 0x3e6e77db,
      0xaed16a4a, 0xd9d65adc, 0x40df0b66, 0x37d83bf0,
      0xa9bcae53, 0xdebb9ec5, 0x47b2cf7f, 0x30b5ffe9,
      0xbdbdf21c, 0xcabac28a, 0x53b39330, 0x24b4a3a6,
      0xbad03605, 0xcdd70693, 0x54de5729, 0x23d967bf,
      0xb3667a2e, 0xc4614ab8, 0x5d681b02, 0x2a6f2b94,
      0xb40bbe37, 0xc30c8ea1, 0x5a05df1b, 0x2d02ef8d
   }; /*crc32Tbl_checkSum*/

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
){
   unsigned long ulByte = 0;
   unsigned short oldUS = 0;

   /*the << 16 in the adler32 is to convert the tow 16 bit
   `  numbers to a single 32 bit number. In this case I am
   `  assining the numbers to their correct little endin
   `  possitions
   */
   unsigned short *oneUS =
      (unsigned short *) &adlerAryUC[2];
   unsigned short *twoUS =
      (unsigned short *) &adlerAryUC[0];

   for(ulByte = 0; ulByte < lenUL; ++ulByte)
   { /*Loop: find adler32 checksum*/
      oldUS = *oneUS;
      *oneUS += dataAryUC[ulByte];

      /*this cases will be rare (about 1 time in in 253)
      `   so branch prediction will skip this most of the
      `   time
      */
      if(*oneUS >= def_adler32_checkSum)
         *oneUS -= def_adler32_checkSum;
      else if(*oneUS < oldUS)
         *oneUS += def_adler32ToUSMax_checkSum;

      oldUS = *twoUS;
      *twoUS += *oneUS;

      /*this cases will be rare (about 1 time in in 253)
      `   so branch prediction will skip this most of the
      `   time
      */
      if(*twoUS >= def_adler32_checkSum)
         *twoUS -= def_adler32_checkSum;
      else if(*twoUS < oldUS)
         *twoUS += def_adler32ToUSMax_checkSum;
   } /*Loop: find adler32 checksum*/

   *oneUS = usToLittle_endin(*oneUS);
   *twoUS = usToLittle_endin(*twoUS);
} /*usAdler32_checkSum*/

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
){
   unsigned short oldUS = 0;

   /*the << 16 in the adler32 is to convert the tow 16 bit
   `  numbers to a single 32 bit number. In this case I am
   `  assining the numbers to their correct little endin
   `  possitions
   */
   unsigned short *oneUS =
      (unsigned short *) &adlerAryUC[2];
   unsigned short *twoUS =
      (unsigned short *) &adlerAryUC[0];

   oldUS = *oneUS;
   *oneUS += byteUC;

   /*this cases will be rare (about 1 time in in 253)
   `   so branch prediction will skip this most of the
   `   time
   */
   if(*oneUS >= def_adler32_checkSum)
      *oneUS -= def_adler32_checkSum;
   else if(*oneUS < oldUS)
      *oneUS += def_adler32ToUSMax_checkSum;

   oldUS = *twoUS;
   *twoUS += *oneUS;

   /*this cases will be rare (about 1 time in in 253)
   `   so branch prediction will skip this most of the
   `   time
   */
   if(*twoUS >= def_adler32_checkSum)
      *twoUS -= def_adler32_checkSum;
   else if(*twoUS < oldUS)
      *twoUS += def_adler32ToUSMax_checkSum;
} /*usAdler32Byte_checkSum*/

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
){
    unsigned short *oneUS = (unsigned short *) &aryUC[2];
    unsigned short *twoUS = (unsigned short *) &aryUC[0];

    *oneUS  = usToLittle_endin(*oneUS);
    *twoUS = usToLittle_endin(*twoUS);
} /*usAdler32Finish_checkSum*/

/*-------------------------------------------------------\
| Fun04: adler32_checkSum
|   - finds adler32 checksum with ints; 32 bit
| Input:
|   - byteUC:
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
){
   unsigned long ulByte = 0;
   signed int oneSI = 0;
   signed int twoSI = 0;

   for(ulByte = 0; ulByte < lenUL; ++ulByte)
   { /*Loop: find adler32 checksum*/
      oneSI += dataAryUC[ulByte];
      oneSI %= def_adler32_checkSum;

      twoSI += oneSI;
      twoSI %= def_adler32_checkSum;
   } /*Loop: find adler32 checksum*/

   return (twoSI << 16) | oneSI;
} /*adler32_checkSum*/

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
){
   *oneSIPtr += byteUC;
   *oneSIPtr %= def_adler32_checkSum;

   *twoSIPtr += *oneSIPtr;
   *twoSIPtr %= def_adler32_checkSum;
} /*adler32Byte_checkSum*/

/*-------------------------------------------------------\
| Fun07: crc32Mac_checkSum
|   - find new crc value after adding next byte
|   - copied from misc0110's libattpng repository, I
|     really do not know whats going on here
| Input:
|   - byte:
|     o byte to convert (
|   - crc32:
|     o current value of input crc
|     o unsigned int is always 32 bit on 32/64 bit OSs
|       or at least I have yet to see a case it is
|       not
| Output:
|   - Returns:
|     o updated crc value
\-------------------------------------------------------*/
/*unsigned int
crc32Mac_checkSum(
   unsigned char byte,
   unsigned int crc32
){
   unsigned int indexUI = crc32;
   indexUI ^= byte;
   indexUI &= 0xff;
   indexUI = crc32Tbl_checkSum[indexUI];
   crc32 >>= 8;
   crc32 ^= indexUI;
   return crc32;
}*/
#define crc32Mac_checkSum(byte, crc32) ( ( crc32Tbl_checkSum[ ((unsigned int) crc32 ^ (unsigned char) (byte)) & 0xff] ^ ( (unsigned int) (crc32) >> 8 ) ))
/*not sure what is going on here, this is the best
`   I can understand it by, but it shows my ingorence
`   - invert: crc32 ^ byte:
`   - to8bit: invert & 255:
`     o converts inverted value to value in lookup table
`   - newVal: crc32Mac_checkSumUL[to8bit]:
`     o find the new crc32 value in the lookup table
`   - clear: crc32 >> 8:
`     o clears old 8 bit value from crc32
`   - newVal ^ clear:
`     o inverts int lookup by remaing crc32 value
*/

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
){
   unsigned int crc32UI = (unsigned int) -1;

   while(lenUL--)
   { /*Loop: get crc32*/
      crc32UI = crc32Mac_checkSum(*aryUC, crc32UI);
      ++aryUC;
   } /*Loop: get crc32*/

   return crc32Finish_checkSum(crc32UI);
} /*crc32_checkSum*/

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
){
   return crc32Mac_checkSum(byteUC, crc32UI);
} /*crc32Byte_checkSum*/

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
){
   unsigned int crc32UI = 0;
   unsigned int uiElm = 0;
   unsigned int uiBit = 0;

   for(uiElm = 0; uiElm < 256; ++uiElm)
   { /*Loop: build crc32 table*/
      crc32UI = uiElm;

      for(uiBit = 0; uiBit < 8; ++uiBit)
      { /*Loop: find crc32 value for uiElm*/
         if(crc32UI & 1)
            crc32UI = (crc32UI >> 1) ^ 0xEDB88320;
         else
             crc32UI >>= 1;
      } /*Loop: find crc32 value for uiElm*/

      crc32AryUI[uiElm] = crc32UI;
   } /*Loop: build crc32 table*/
} /*mkCrc32Tbl_checkSum*/

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
