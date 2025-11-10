/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' kmerBit SOF: Start Of File
'   - has functions to converts a nucleotide index
'     sequence (seqToIndex_alnSet) into to a kmer sequence
'   o header:
'     - included libraries
'   o fun01: seqToKmer_kmerFun
'     - convert sequence to a kmer list
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

#include "kmerFun.h"
#include "kmerBit.h" /*has tables/macros for conversion*/

/*-------------------------------------------------------\
| Fun01: seqToKmer_kmerFun
|   - convert sequence to a kmer list
| Input:
|   - seqStr:
|     o c-string with sequence to convert
|   - lenSeqSI:
|     o length of sequence to convert
|   - kmerArySIPtr:
|     o signed int array pionter to hold kmer sequence
|   - sizeSIPtr:
|     o signed int pointer with current kmerArySIPtr size
|   - kmerLenUC:
|     o length of one kmer (index 1)
| Output:
|   - Modifies:
|     o kmerArySIPtr to have forward kmers (resized if
|       needed)
|     o sizeSIPtr to have resized array length
|   - Returns:
|     o number kmers in sequence
|     o 0 for memory error
\-------------------------------------------------------*/
signed int
seqToKmer_kmerFun(
   signed char *seqStr,        /*sequence to convert*/
   signed int lenSeqSI,        /*sequence length*/
   signed int **kmerArySIPtr,  /*gets kmer sequence*/
   signed int *sizeSIPtr,      /*current array sizes*/
   unsigned char kmerLenUC     /*size of one kmer*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - get forward and reverse kmer sequences
   '   o fun01 sec01:
   '     - variable declarations
   '   o fun01 sec02:
   '     - allocate memory for arrays
   '   o fun01 sec03:
   '     - add kmers to arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned char ntUC = 0;
   signed int siNt = 0;

   unsigned long kmerUL = 0; /*forward kmer*/
   unsigned long lenUL = 0;  /*number kmers*/
   unsigned long maskUL = mkMask_kmerBit(kmerLenUC);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - allocate memory for arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(*sizeSIPtr < lenSeqSI)
   { /*If: need more memory*/
      if(*kmerArySIPtr)
         free(*kmerArySIPtr);
      *kmerArySIPtr = 0;
   } /*If: need more memory*/

   if(! *kmerArySIPtr)
   { /*If: need to allocate memory*/
      *kmerArySIPtr =
         malloc( (lenSeqSI + 8) * sizeof(signed int) );
      if(! *kmerArySIPtr)
         goto memErr_fun01_sec04;
   } /*If: need to allocate memory*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - add kmers to array
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   --kmerLenUC; /*convert to index 0*/

   for(
      siNt = 0;
      siNt < lenSeqSI;
      ++siNt
   ){ /*Loop: add kmers*/

      ntUC =
         alnNtTo_kmerBit[ (unsigned char) seqStr[siNt] ];

      kmerUL =
         ntBitToKmer_kmerBit(
            ntUC,
            kmerUL,
            maskUL
         );

      if(ntUC < def_anonNt_kmerBit)
      { /*If: no errors/anonymous bases*/
         ++lenUL; /*kmer length*/

         if(siNt < kmerLenUC)
            ; /*on first few bases*/
         else if(lenUL < kmerLenUC)
            (*kmerArySIPtr)[siNt - kmerLenUC] =
                def_noKmer_kmerBit;
            /*had anonymous base*/
         else
            (*kmerArySIPtr)[siNt - kmerLenUC] =
               (signed int) kmerUL;
      } /*If: no errors/anonymous bases*/

      else
         lenUL = 0;
   }  /*Loop: add kmers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return lenSeqSI - kmerLenUC;

   memErr_fun01_sec04:;
      return 0;
} /*seqToKmer_kmerFun*/

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
