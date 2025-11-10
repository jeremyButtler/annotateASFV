/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' kmerBit SOF: Start Of File
'   - has table and macros to converts a set of three
'     nucleotide index's (seqToIndex_alnSet) to a kmer
'   o .h tbl01: alnNtTo_kmerBit
'     - converts an nucleotide alignment code form
'       alnSetStruct.h to an two bit value, with an 3rd
'       bit being used for anonymous bases and the 4th bit
'       for errors
'   o .h fun01: mkMask_kmerBit
'     - makes a mask for clearing extra bits from a kmer
'   o .h fun02: ntIndexToKmer_kmerBit
'     - adds a nucleotide (in index format
'      [seqToIndex_alnSet]) to a kmer
'   o .h fun03: ntBitToKmer_kmerBit
'     - adds a nucleotide bit (by alNtTo_kmerBit) to kmer
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifndef KMER_BIT_H
#define KMER_BIT_H

/*if you change these, then also change these in
`  kmerFun.h
*/
#define def_bitsPerKmer_kmerBit 2 /*do not change*/
#define def_endKmers_kmerBit -2
#define def_noKmer_kmerBit -1

#define def_a_kmerBit 0
#define def_c_kmerBit 1
#define def_g_kmerBit 2
#define def_t_kmerBit 3
#define def_anonNt_kmerBit 4 /*anonmyous nucleotide*/
#define def_noNt_kmerBit 8

/*-------------------------------------------------------\
| Tbl01: alnNtTo_kmerBit
|   - converts an nucleotide alignment code form
|     alnSetStruct.h to an two bit value, with an 3rd
|     bit being used for anonymous bases and the 4th bit
|     for errors
|   - this will generate a warning if not used. Bit of
|     a pain, but need the table to be public
\-------------------------------------------------------*/
static unsigned char alnNtTo_kmerBit[] =
{
   def_noNt_kmerBit,    /*00 = 64 = ,*/  
   def_a_kmerBit,       /*01 = 65 = A/a*/  
   def_anonNt_kmerBit,  /*02 = 66 = B/b*/  
   def_c_kmerBit,       /*03 = 67 = C/c*/  
   def_anonNt_kmerBit,  /*04 = 68 = D/d*/  
   def_noNt_kmerBit,    /*05 = 69 = E/e*/  
   def_noNt_kmerBit,    /*06 = 70 = F/f*/  
   def_g_kmerBit,       /*07 = 71 = G/g*/  
   def_anonNt_kmerBit,  /*08 = 72 = H/h*/  
   def_noNt_kmerBit,    /*09 = 73 = I/i*/  
   def_noNt_kmerBit,    /*10 = 74 = J/j*/  
   def_anonNt_kmerBit,  /*11 = 75 = K/k*/  
   def_noNt_kmerBit,    /*12 = 76 = L/l*/  
   def_anonNt_kmerBit,  /*13 = 77 = M/m*/  
   def_anonNt_kmerBit,  /*14 = 78 = N/n*/  
   def_noNt_kmerBit,    /*15 = 79 = O/o*/  
   def_noNt_kmerBit,    /*16 = 80 = P/p*/  
   def_noNt_kmerBit,    /*17 = 81 = Q/q*/  
   def_anonNt_kmerBit,  /*18 = 82 = R/r*/  
   def_anonNt_kmerBit,  /*19 = 83 = S/s*/  
   def_t_kmerBit,       /*20 = 84 = T/t*/  
   def_t_kmerBit,       /*21 = 85 = U/u*/  
   def_anonNt_kmerBit,  /*22 = 86 = V/v*/  
   def_anonNt_kmerBit,  /*23 = 87 = W/w*/  
   def_anonNt_kmerBit,  /*24 = 88 = X/x (amino acids)*/  
   def_anonNt_kmerBit,  /*25 = 89 = Y/y*/  
   def_noNt_kmerBit     /*26 = 90 = Z/z*/  
}; /*alnNtTo_kmerBit*/

/*-------------------------------------------------------\
| Fun01: mkMask_kmerBit
|   - makes a mask for clearing extra bits from a kmer
| Input:
|   - lenKmerMac:
|     o number nucleotides in one kmer
| Output:
|   - Returns:
|     o mask to clear kmer (0's for unused bits and 1's
|       for used bits)
\-------------------------------------------------------*/
#define mkMask_kmerBit(lenKmerMac) ( ((unsigned long) -1) >> (( sizeof(unsigned long) << 3 ) - ( (lenKmerMac) * def_bitsPerKmer_kmerBit )) )
/*Logic:
`  - kmerBits: lenKmerMac * def_bitsPerKmer_kmerBit:
`    o gets number of bits needed to fill one kmer
`  - sizeLong: sizeof(unsigned long) << 3
`    o number of bits in unsigned long
`  - emptyBits: sizeLong - kmerBits:
`    o number of extra bits in unsigned long (not used)
`  - ( (unsinged long) -1 ) >> emptyBits
`    o shifts -1 till all unsed bits are set to 0, and all
`      used bits are set to 1
`     
*/

/*-------------------------------------------------------\
| Fun02: ntIndexToKmer_kmerBit
|   - adds a nucleotide (in index format
|    [seqToIndex_alnSet]) to a kmer
| Input:
|   - ntMac:
|     o nucleotide to add
|   - kmerMac:
|     o kmer with previous nucleotides
|   - maskMac:
|     o mask to clear extr bits with (use mkMask_kmerBit)
| Output:
|   - Returns:
|     o kmer with added bits (does not modify kmerMac)
\-------------------------------------------------------*/
#define ntIndexToKmer_kmerBit(ntMac, kmerMac, maskMac) ( ( ((kmerMac) << def_bitsPerKmer_kmerBit) | (alnNtTo_kmerBit[ (unsigned char) (ntMac) ]) ) & (maskMac) )
/*Logic:
`   - mkRoom: kmerMac << def_bitsPerKmer_kmerBit:
`     o add room for new bit
`   - getBit: alnNtTo_kmerBit[(unsigned char) ntMac]:
`     o convert nucleotide index to bits for kmer
`   - kmer: mkRoom | getBit:
`     o add converted nucleotide index to kmer
`   - kmer & maskMac:
`     o clear any extra bits (not in kmer)
*/

/*-------------------------------------------------------\
| Fun03: ntBitToKmer_kmerBit
|   - adds a nucleotide bit (by alNtTo_kmerBit) to a kmer
| Input:
|   - ntBitMac:
|     o nucleotide bit to add; covert by alnNtTo_kmerBit
|   - kmerMac:
|     o kmer with previous nucleotides
|   - maskMac:
|     o mask to clear extr bits with (use mkMask_kmerBit)
| Output:
|   - Returns:
|     o kmer with added bits (does not modify kmerMac)
\-------------------------------------------------------*/
#define ntBitToKmer_kmerBit(ntBitMac, kmerMac, maskMac) ( ( ((kmerMac) << def_bitsPerKmer_kmerBit) | ((unsigned char) (ntBitMac)) ) & (maskMac) )
/*Logic:
`   - mkRoom: kmerMac << def_bitsPerKmer_kmerBit:
`     o add room for new bit
`   - kmer: mkRoom | ntBitMac:
`     o add converted nucleotide index to kmer
`   - kmer & maskMac:
`     o clear any extra bits (not in kmer)
*/

#endif

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
