/*#########################################################
# Name: vecWrapAvx512Fun.h
# Use:
#  - Holds macro wrappers and variables for AVX512
#    intrinsics
# Note:
#  - Do not use the min/max for unsigned 32 ints if you
#    are planning for SSE2 support. At best (no load costs)
#    it has the same efficany as doing the equivlent number
#    of min/max's.
#########################################################*/

#ifndef VECTWRAPAVX512_H
#define VECTWRAPAVX512_H

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' vectWrapAvx512Fun.h TOC: Table Of Contents
'  - Holds macro wrapper for AVX512 intrinsics
'  o sec-01:
'    - Definitions and variable declerations
'  o sec-02:
'    - Macros wrapping each intrinsic
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-01:
^  - Definitions and variable declerations
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#include <immintrin.h>
#include <stdint.h>

#define defVectBits 512
#define defVectBytes 64
#define defNum8BitElms 64
#define defNum16BitElms 32
#define defNum32BitElms 16
#define defMvMask8 0xffffffffffffffff
#define defMvMask16 0xffffffff
#define defMvMask32 0xffff

typedef __m512i vectI8;  /*vector of bytes (8bits)*/
typedef __m512i vectI16; /*vector of shorts (16 bits)*/
typedef __m512i vectI32; /*vector of ints (32 bits)*/

typedef __m512i vectU8;  /*Prefer unsigned bytes*/
typedef __m512i vectU16; /*Prefer unsigned shorts*/
typedef __m512i vectU32; /*Prefer unsinged ints*/

typedef __mask64 mask8;  /*mask of 8 bit values*/
typedef __mask32 mask16; /*Mask of 16 bit values*/
typedef __mask16 mask32; /*mask of 32 bit values*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-02:
^  - Macros
^  o sec-02 sub-01:
^    - Input/output functions (load, set, store)
^  o sec-02 sub-02:
^    - Output functions
^  o sec-02 sub-03:
^    - Comparison functions (if)
^  o sec-02 sub-04:
^    - Logical functions (andNot, and, or, xor)
^  o sec-02 sub-05:
^    - Bit manipulation functions (shifts)
^  o sec-02 sub-06:
^    - Math functions
^  o sec-02 sub-07:
^    - Casts [for NEON & masks] and Conversions
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*********************************************************\
* Sec-02 Sub-01:
*  - Input functions
*  o sec-02 sub-01 cat-01:
*    - aligned loading
*  o sec-02 sub-01 cat-02:
*    - Unaligned loading
*  o sec-02 sub-01 cat-03:
*    - make zero vectors
*  o sec-02 sub-01 cat-04:
*    - Make vectors of one element
*  o sec-02 sub-01 cat-05:
*    - Insert an element into an vector
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-01:
+  - Aligned loading
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define load_I8Ary_retVectI8(arrayI8) \
  (_mm512_load_si512((__m512i *) (arrayI8));)

#define load_I16Ary_retVectI16(arrayI16) \
  (_mm512_load_si512((__m512i *) (arrayI16));)

#define load_I32Ary_retVectI32(arrayI32) \
  (_mm512_load_si512((__m512i *) (arrayI32));)


#define load_U8Ary_retVectU8(arrayU8) \
  (_mm512_load_si512((__m512i *) (arrayU8));)

#define load_U16Ary_retVectU16(arrayU16) \
  (_mm512_load_si512((__m512i *) (arrayU16));)

#define load_U32Ary_retVectU32(arrayU32) \
  (_mm512_load_si512((__m512i *) (arrayU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-02:
+  - Unaligned loading
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define loadu_I8Ary_retVectI8(arrayI8) \
  (_mm512_loadu_si512((__m512i *) (arrayI8));)

#define loadu_I16Ary_retVectI16(arrayI16) \
  (_mm512_loadu_si512((__m512i *) (arrayI16));)

#define loadu_I32Ary_retVectI32(arrayI32) \
  (_mm512_loadu_si512((__m512i *) (arrayI32));)


#define loadu_U8Ary_retVectU8(arrayU8) \
  (_mm512_loadu_si512((__m512i *) (arrayU8));)

#define loadu_U16Ary_retVectU16(arrayU16) \
  (_mm512_loadu_si512((__m512i *) (arrayU16));)

#define loadu_U32Ary_retVectU32(arrayU32) \
  (_mm512_loadu_si512((__m512i *) (arrayU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-03:
+  - make zero vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define zero_retVectI8() (_mm512_setzero_si512())
#define zero_retVectI16()(_mm512_setzero_si512())
#define zero_retVectI32()(_mm512_setzero_si512())

#define zero_retVectU8() (_mm512_setzero_si512())
#define zero_retVectU16()(_mm512_setzero_si512())
#define zero_retVectU32()(_mm512_setzero_si512())

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-04:
+  - Make vectors of one element
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define set1_I8_retVectI8(valC) (mm512_set1_epi8((valC)))

#define set1_I16_retVectI16(valI16) \
   (mm512_set1_epi16((valI16)))

#define set1_I32_retVectI32(valI32) \
   (_mm512_set1_epi32((valI32)))


#define set1_U8_retVectU8(valUC)(mm512_set1_epi8((valUC)))

#define set1_U16_retVectU16(valC) \
   (mm512_set1_epi16((valU16)))

#define set1_U32_retVectU32(valU32) \
   (_mm512_set1_epi32((valU32)))


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-05:
+  - Insert an element into an vector
+  - https://stackoverflow.com/questions/58303958/how-to-implement-16-and-32-bit-integer-insert-and-extract-operations-with-avx-512
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*#define insert_I8_retVectI8(insI8, posI){\
   _mm512_mask_set1_epi8(insVectI8,1UL<<posI,insI8);\
}*/ /*No SSE2 support (no easy conversion either)*/

#define insert_I16_retVectI16(insVectI16, insI16, posI)\
(   _mm512_mask_set1_epi16(\
      (insVectI16),\
      1UL<<(posI),\
      (insI16));\
)

#define insert_I32_retVectI32(insVectI32, insI32, posI)\
(   _mm512_mask_set1_epi32(\
      (insVectI32),\
      1UL<<(posI),\
      (insI32));\
)


/*#define insert_U8_retVectU8(insU8, posU){\
   _mm512_mask_set1_epi8(insVectU8,1UL<<posU,insU8);\
}*/ /*No SSE2 support (no easy conversion either)*/


#define insert_U16_retVectU16(insVectU16, insU16, posU)\
(   _mm512_mask_set1_epi16(\
      (insVectU16),\
      1UL << (posU),\
      (insU16));\
)

#define insert_U32_retVectU32(insVectU32, insU32, posU)\
(   _mm512_mask_set1_epi32(\
       (insVectU32),\
       1UL << (posU),\
       (insU32)); \
)

/*********************************************************\
* Sec-02 Sub-02:
*  - Output functions
*  o sec-02 sub-02 cat-01:
*    - Store vector ouput into an array (aligned)
*  o sec-02 sub-02 cat-02:
*    - Store vector ouput into an array (unaligned)
*  o sec-02 sub-02 cat-04:
*    - Store masks into longs
*  o sec-02 sub-02 cat-04:
*    - Extract an element
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-01:
+  - Store vector ouput into an array (aligned)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define store_vectI8_retAryI8(retAry, inVectI8) \
  (_mm512_store_si512((__m512i *) (retAry), (inVectI8));)

#define store_vectI16_retAryI16(retAry, inVectI16) \
  (_mm512_store_si512((__m512i *) (retAry), (inVectI16));)

#define store_vectI32_retAryI32(retAry, inVectI32) \
  (_mm512_store_si512((__m512i *) (retAry), (inVectI32));)


#define store_vectU8_retAryU8(retAry, inVectU8) \
  (_mm512_store_si512((__m512i *) (retAry), (inVectU8));)

#define store_vectU16_retAryU16(retAry, inVectU16) \
  (_mm512_store_si512((__m512i *) (retAry), (inVectU16));)

#define store_vectU32_retAryU32(retAry, inVectU32) \
  (_mm512_store_si512((__m512i *) (retAry), (inVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-02:
+  - Store vector ouput into an array (unaligned)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define storeu_vectI8_retAryI8(retAry, inVectI8) \
  (_mm512_storeu_si512((__m512i *) (retAry), (inVectI8));)

#define storeu_vectI16_retAryI16(retAry, inVectI16) \
  (_mm512_storeu_si512((__m512i *) (retAry), (inVectI16));)

#define storeu_vectI32_retAryI32(retAry, inVectI32) \
  (_mm512_storeu_si512((__m512i *) (retAry), (inVectI32));)


#define storeu_vectU8_retAryU8(retAry, inVectU8) \
  (_mm512_storeu_si512((__m512i *) (retAry), (inVectU8));)

#define storeu_vectU16_retAryU16(retAry, inVectU16) \
  (_mm512_storeu_si512((__m512i *) (retAry), (inVectU16));)

#define storeu_vectU32_retAryU32(retAry, inVectU32) \
  (_mm512_storeu_si512((__m512i *) (retAry), (inVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-03:
+  - Store masks into longs
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define store_mask8_retUL(retUL, mask8) \
  (_store_mask64((__mask64 *) (retUL), (mask8));)

#define store_mask16_retUL(retUL, mask16) \
  (_store_mask32((__mask32 *) (retUL), (mask16));)

#define store_mask32_retUL(retUL, mask32) \
  (_store_mask16((__mask16 *) (retUL), (mask32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-04:
+  - Extract an element
+  - This code idea was modified from
+  - https://stackoverflow.com/questions/75997936/efficiently-extract-single-double-element-from-avx-512-vector
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*Not very efficent, but might work*/
#define extract_vectI8_I8(inVectI8, charToExtS)({\
   __mm512i shuffleVI;\
   \
   shuffleVI = \
      _mm512_casti128_si512(/*Converts _mm128i to_mm512i*/\
         _mm_cvtsi32_si128((charToExtS) >> 2)\
         /*This puts an int to the lowest position in an
         ` _mm128i vector (4 chars per int) Latency 1 to 2
         */\
   )\
   \
   _mm512_cvtsi512_si32(\
         /*Extracts lowest integer (latency 2 to 3)*/\
      _mm512_permutexvar_epi8(shuffleVI, (inVectI8))\
        /*Moves target to lowest position (laten 4 to 6)*/\
   ) /*Extract the lowest 32 bit value*/\
   & 0xFF << (( ( (charToExtS) >> 4 ) - 1 ) << 3)\
   /*0xFF is the first 8 bits
   `    charToExtS >> 4
   `       Gives the position of the char in the int
   `    (charToExtS >> 4) - 1
   `       Converts the chars position to index 0
   `    ((charToExtS >> 4) - 1) << 3
   `       Converts the position to 0, 8, 16, or 24
   `    0xFF (((charToExtS >> 4) - 1) << 3)
   `       shifts 0xFF by 0, 8, 16, or 24 bits. This
   `       is the index of my target char in an integer.
   `       All I need to do is clear the nontaget bits
   `       with an & operation.
   */\
})

#define extract_vectU8_U8(inVectU8, charToExtS)({\
   __mm512i shuffleVI;\
   \
   shuffleVI = \
      _mm512_casti128_si512(/*Converts _mm128i to_mm512i*/\
         _mm_cvtsi32_si128((charToExtS) >> 2)\
         /*This puts an int to the lowest position in an
         ` _mm128i vector (4 chars per int)
         */\
   )\
   \
   _mm512_cvtsi512_si32( /*Extracts lowest integer*/\
      _mm512_permutexvar_epi8(shuffleVI, (inVectU8))\
      /*Moves target to the lowest position*/\
   ) /*Extract the lowest 32 bit value*/\
   & 0xFF << (( ( (charToExtS) >> 4 ) - 1 ) << 3)\
   /*0xFF is the first 8 bits
   `    charToExtS >> 4
   `       Gives the position of the char in the int
   `    (charToExtS >> 4) - 1
   `       Converts the chars position to index 0
   `    ((charToExtS >> 4) - 1) << 3
   `       Converts the position to 0, 8, 16, or 24
   `    0xFF (((charToExtS >> 4) - 1) << 3)
   `       shifts 0xFF by 0, 8, 16, or 24 bits. This
   `       is the index of my target char in an integer.
   `       All I need to do is clear the nontaget bits
   `       with an & operation.
   */\
})


#define extract_vectI16_I16(inVectI16, shortToExtS)({\
   __mm512i shuffleVI;\
   \
   shuffleVI = \
      _mm512_casti128_si512(/*Converts _mm128i to_mm512i*/\
         _mm_cvtsi32_si128((shortToExtS) >> 1)\
         /*This puts an int to the lowest position in an
         ` _mm128i vector (2 shorts per int)
         ` The >> 1 converts the shorts position to int
         */\
   )\
   \
   (int16_t) _mm512_cvtsi512_si32(\
         /*Extracts lowest integer*/\
      _mm512_permutexvar_epi16(shuffleVI, (inVectI16))\
      /*Moves target to the lowest position*/\
   ) /*Extract the lowest 32 bit value*/\
   >> (((shortToExtS) & 1) >> 4)\
   /*
   `  shortToExtS & 1 tells me if this was odd or even
   `  (shortToExtS & 1) >> 4
   `      is 16 if shortToExtS is odd, 0 if even
   `   extracted number >> 16 is the last short (odd),
   `   (int16_t) extracted number >> 0 clears the frist
   `     bits
   `   I used to do an 0xFFFFFFFF << 16 or 0, but this
   `     is a bit cleaner
   */\
})

#define extract_vectU16_U16(inVectU16, shortToExtS)({\
   __mm512i shuffleVI;\
   \
   shuffleVI = \
      _mm512_casti128_si512(/*Converts _mm128i to_mm512i*/\
         _mm_cvtsi32_si128((shortToExtS) >> 1)\
         /*This puts an int to the lowest position in an
         ` _mm128i vector (2 shorts per int)
         ` The >> 1 converts the shorts position to int
         */\
   )\
   \
   (uint16_t) _mm512_cvtsi512_si32(\
          /*Extracts lowest integer*/\
      _mm512_permutexvar_epi16(shuffleVI, (inVectU16))\
      /*Moves target to the lowest position*/\
   ) /*Extract the lowest 32 bit value*/\
   >> (((shortToExtS) & 1) >> 4)\
   /*
   `  shortToExtS & 1 tells me if this was odd or even
   `  (shortToExtS & 1) >> 4
   `      is 16 if shortToExtS is odd, 0 if even
   `   extracted number >> 16 is the last short (odd),
   `   (int16_t) extracted number >> 0 clears the frist
   `     bits
   `   I used to do an 0xFFFFFFFF << 16 or 0, but this
   `     is a bit cleaner
   */\
})


#define extract_vectI32_I32(inVectI32, intToExt)({\
   __mm512i shuffleVI;\
   \
   shuffleVI = \
      _mm512_casti128_si512(/*Converts _mm128i to_mm512i*/\
         _mm_cvtsi32_si128((intToExt))\
         /*This puts intToExt to the lowest position in an
         ` __mm128i vector
         */\
   )/*Convert intToExt to a vector; Casting has no cost*/\
   \
   _mm512_cvtsi512_si32( /*Extracts lowest integer*/\
      _mm512_permutexvar_epi32(shuffleVI, (inVectI32))\
         /*Moves target to the lowest position*/\
   ) /*Extract the lowest 32 bit value*/\
})


#define extract_vectU32_U32(inVectU32, intToExt)({\
   __mm512i shuffleVI;\
   \
   shuffleVI = \
      _mm512_casti128_si512(/*Converts _mm128i to_mm512i*/\
         _mm_cvtsi32_si128((intToExt))\
         /*This puts intToExt to the lowest position in an
         ` __mm128i vector
         */\
   )/*Convert intToExt to a vector; Casting has no cost*/\
   \
   _mm512_cvtsi512_si32( /*Extracts lowest integer*/\
      _mm512_permutexvar_epi32(shuffleVI, (inVectU32))\
         /*Moves target to the lowest position*/\
   ) /*Extract the lowest 32 bit value*/\
})

/*********************************************************\
* Sec-02 Sub-03:
*  - Comparison functions (if)
*  o sec-02 sub-03 cat-01:
*    - Comparisons (return masks)
*  o sec-02 sub-03 cat-02:
*    - Unsinged equal comparisions (returns masks)
*  o sec-02 sub-03 cat-03:
*    - Unsinged less than comparisions (less efficent;mask)
*  o sec-02 sub-03 cat-04:
*    - Unsinged greater than comparisions (less efficent)
*  o sec-02 sub-03 cat-05:
*    - Fix differences in population counts (total 1's)
*    - This is here because I am also supporting SSE2 and
*      NEON.
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-01:
+  - Comparisons (return masks)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*8-bit comparisions; returns a mask8*/
#define cmpeq_vectI8_retMask8(xVectI8, yVectI8) \
  ((mask8) _mm512_cmpeq_epi8_mask((xVectI8),(yVectI8));)

#define cmpgt_vectI8_retMask8(xVectI8, yVectI8) \
  ((mask8) _mm512_cmpgt_epi8_mask((xVectI8),(yVectI8));)

#define cmplt_vectI8_retMask8(xVectI8, yVectI8) \
  ((mask8) _mm512_cmplt_epi8_mask((xVectI8),(yVectI8));)

/*16-bit comparisions; returns a mask16*/

#define cmpeq_vectI16_retMask16(xVectI16, yVectI16) \
 ((mask16) _mm512_cmpeq_epi16_mask((xVectI16),(yVectI16));)

#define cmpgt_vectI16_retMask16(xVectI16, yVectI16) \
 ((mask16) _mm512_cmpgt_epi16_mask((xVectI16),(yVectI16));)

#define cmplt_vectI16_retMask16(xVectI16, yVectI16) \
 ((mask16) _mm512_cmplt_epi16_mask((xVectI16),(yVectI16));)

/*32-bit comparisions; returns a mask32*/

#define cmpeq_vectI32_retMask32(xVectI32, yVectI32) \
 ((mask32) _mm512_cmpeq_epi32_mask((xVectI32),(yVectI32));)

#define cmpgt_vectI32_retMask32(xVectI32, yVectI32) \
 ((mask32) _mm512_cmpgt_epi32_mask((xVectI32),(yVectI32));)

#define cmplt_vectI32_retMask32(xVectI32, yVectI32) \
 ((mask32) _mm512_cmplt_epi32_mask((xVectI32),(yVectI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-02:
+  - Unsinged equal comparisions (returns masks)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cmpeq_vectU8_retMask8(xVectI8, yVectI8) \
   ((mask8) _mm512_cmpeq_epu8_mask((xVectI8),(yVectI8));)

#define cmpeq_vectU16_retMask16(xVectI16, yVectI16) \
 ((mask16) _mm512_cmpeq_epu16_mask((xVectI16),(yVectI16));)

#define cmpeq_vectU32_retMask32(xVectI32, yVectI32) \
 ((mask32) _mm512_cmpeq_epu32_mask((xVectI32),(yVectI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-03:
+  - Unsinged < comparisions (non-AVX512 less efficent)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cmplt_vectU8_retMask8(xVectI8, yVectI8) \
   ((mask8) _mm512_cmplt_epu8_mask((xVectI8),(yVectI8));)

#define cmplt_vectU16_retMask16(xVectI16, yVectI16) \
 ((mask16) _mm512_cmplt_epu16_mask((xVectI16),(yVectI16));)

#define cmplt_vectU32_retMask32(xVectI32, yVectI32) \
 ((mask32) _mm512_cmplt_epu32_mask((xVectI32),(yVectI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-04:
+  - Unsinged > comparisions (non-AVX512 less efficent)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cmpgt_vectU8_retMask8(xVectI8, yVectI8) \
   ((mask8) _mm512_cmpgt_epu8_mask((xVectI8),(yVectI8));)

#define cmpgt_vectU16_retMask16(xVectI16, yVectI16) \
 ((mask16) _mm512_cmpgt_epu16_mask((xVectI16),(yVectI16));)

#define cmpgt_vectU32_retMask32(xVectI32, yVectI32) \
 ((mask32) _mm512_cmpgt_epu32_mask((xVectI32),(yVectI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-02:
+  - Fix differences in population counts (total 1's)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* This is here to support SSE2, AVX2, and NEON. 
`  The returned mask in AVX512 has only one bit per
`  data type in comparison
*/
#define fix_mask8_popcount(inUL) ((inUL))
#define fix_mask16_popcount(inUL) ((inUL))
#define fix_mask32_popcount(inUL) ((inUL))

/*********************************************************\
* Sec-02 Sub-04:
*  - Logical functions (andNot, and, or, xor)
*  o sec-02 sub-02 cat-01:
*    - and not functions
*  o sec-02 sub-02 cat-02:
*    - and functions
*  o sec-02 sub-02 cat-03:
*    - or functions
*  o sec-02 sub-02 cat-04:
*    - xor functions
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-01:
+  - and not functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define andNot_vectI8_retVectI8(notVectI8, inVectI8) \
   (_mm512_andnot_si512((notVectI8), (inVectI8));)

#define andNot_vectU8_retVectU8(notVectU8, inVectU8) \
   (_mm512_andnot_si512((notVectU8), (inVectU8));)

#define andNot_vectI16_retVectI16(notVectI16, inVectI16) \
   (_mm512_andnot_si512((notVectI16), (inVectI16));)

#define andNot_vectU16_retVectU16(notVectU16, inVectU16) \
   (_mm512_andnot_si512((notVectU16), (inVectU16));)

#define andNot_vectI32_retVectI32(notVectI32, inVectI32) \
   (_mm512_andnot_si512((notVectI32), (inVectI32));)

#define andNot_vectU32_retVectU32(notVectU32, inVectU32) \
   (_mm512_andnot_si512((notVectU32), (inVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-02:
+  - and functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define and_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm512_and_si512((xVectI8), (yVectI8));)

#define and_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm512_and_si512((xVectU8), (yVectU8));)

#define and_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm512_and_si512((xVectI16), (yVectI16));)

#define and_vectU16_retVectU16(xVectU16, yVectU16) \
   (_mm512_and_si512((xVectU16), (yVectU16));)

#define and_vectI32_retVectI32(xVectI32, yVectI32) \
   (_mm512_and_si512((xVectI32), (yVectI32));)

#define and_vectU32_retVectU32(xVectU32, yVectU32) \
   (_mm512_and_si512((xVectU32), (yVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-03:
+  - or functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define or_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm512_or_si512((xVectI8), (yVectI8));)

#define or_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm512_or_si512((xVectU8), (yVectU8));)

#define or_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm512_or_si512((xVectI16), (yVectI16));)

#define or_vectU16_retVectU16(xVectU16, yVectU16) \
   (_mm512_or_si512((xVectU16), (yVectU16));)

#define or_vectI32_retVectI32(xVectI32, yVectI32) \
   (_mm512_or_si512((xVectI32), (yVectI32));)

#define or_vectU32_retVectU32(xVectU32, yVectU32) \
   (_mm512_or_si512((xVectU32), (yVectU32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-04:
+  - xor functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define xor_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm512_xor_si512((xVectI8), (yVectI8));)

#define xor_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm512_xor_si512((xVectU8), (yVectU8));)

#define xor_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm512_xor_si512((xVectI16), (yVectI16));)

#define xor_vectU16_retVectU16(xVectU16, yVectU16) \
   (_mm512_xor_si512((xVectU16), (yVectU16));)

#define xor_vectI32_retVectI32(xVectI32, yVectI32) \
   (_mm512_xor_si512((xVectI32), (yVectI32));)

#define xor_vectU32_retVectU32(xVectU32, yVectU32) \
   (_mm512_xor_si512((xVectU32), (yVectU32));)

/*********************************************************\
* Sec-02 Sub-05:
*  - Bit manipulation functions (shifts)
*  o sec-02 sub-05 cat-01:
*    - Shift numbers in vectors right by x bits
*  o sec-02 sub-05 cat-03:
*    - Shift a vector right by x bytes
*  o sec-02 sub-05 cat-02:
*    - Shift numbers in vectors left by x bits
*  o sec-02 sub-05 cat-04:
*    - Shift a vector left by x bytes
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-01:
+  - Shift each element right by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sr_vectI8_retVectI8(inVectI8, numBitsI)\
    (_mm512_srli_epi8((inVectI8), (numBitsI));)

#define sr_vectU8_retVectU8(inVectU8, numBitsI)\
    (_mm512_srli_epi8((inVectU8), (numBitsI));)

#define sr_vectI16_retVectI16(inVectI16, numBitsI)\
    (_mm512_srli_epi16((inVectI16), (numBitsI));)

#define sr_vectU16_retVectU16(inVectU16, numBitsI)\
    (_mm512_srli_epi16((inVectU16), (numBitsI));)

#define sr_vectI32_retVectI32(inVectI32, numBitsI)\
    (_mm512_srli_epi32((inVectI32), (numBitsI));)

#define sr_vectU32_retVectU32(inVectU32, numBitsI)\
    (_mm512_srli_epi32((inVectU32), (numBitsI));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-02:
+  - Shift a vector right by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define srvect_vectI8_retVectI8(inVectI8, numBytesI)\
    (_mm512_srli_si512((inVectI8), (numBytesI));)

#define srvect_vectU8_retVectU8(inVectU8, numBytesI)\
    (_mm512_srli_si512((inVectU8), (numBytesI));)

#define srvect_vectI16_retVectI16(inVectI16, numShortsI)\
    (_mm512_srli_si512((inVectI16), (numShortsI) << 1);)

#define srvect_vectU16_retVectU16(inVectU16, numShortsI)\
    (_mm512_srli_si512((inVectU16), (numShortsI) << 1);)

#define srvect_vectI32_retVectI32(inVectI32, numIntsI)\
    (_mm512_srli_si512((inVectI32), (numIntsI) << 2);)

#define srvect_vectU32_retVectU32(inVectU32, numIntsI)\
    (_mm512_srli_si512((inVectU32), (numIntsI) << 2);)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-03:
+  - Shift numbers in vectors left by x bits
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sl_vectI8_retVectI8(inVectI8, numBitsI)\
    (_mm512_slli_epi8((inVectI8), (numBitsI));)

#define sl_vectU8_retVectU8(inVectU8, numBitsI)\
    (_mm512_slli_epi8((inVectU8), (numBitsI));)

#define sl_vectI16_retVectI16(inVectI16, numBitsI)\
    (_mm512_slli_epi16((inVectI16), (numBitsI));)

#define sl_vectU16_retVectU16(inVectU16, numBitsI)\
    (_mm512_slli_epi16((inVectU16), (numBitsI));)

#define sl_vectI32_retVectI32(inVectI32, numBitsI)\
    (_mm512_slli_epi32((inVectI32), (numBitsI));)

#define sl_vectU32_retVectU32(inVectU32, numBitsI)\
    (_mm512_slli_epi32((inVectU32), (numBitsI));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-04:
+  - Shift a vector left by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define slvect_vectI8_retVectI8(inVectI8, numBytesI)\
    (_mm512_slli_si512((inVectI8), (numBytesI));)

#define slvect_vectU8_retVectU8(inVectU8, numBytesI)\
    (_mm512_slli_si512((inVectU8), (numBytesI));)

#define slvect_vectI16_retVectI16(inVectI16, numShortsI)\
    (_mm512_slli_si512((inVectI16), (numShortsI) << 1);)

#define slvect_vectU16_retVectU16(inVectU16, numShortsI)\
    (_mm512_slli_si512((inVectU16), (numShortsI) << 1);)

#define slvect_vectI32_retVectI32(inVectI32, numIntsI)\
    (_mm512_slli_si512((inVectI32), (numIntsI) << 2);)

#define slvect_vectU32_retVectU32(inVectU32, numIntsI)\
    (_mm512_slli_si512((inVectU32), (numIntsI) << 2);)

/*********************************************************\
* Sec-02 Sub-06:
*  - Math functions
*  o sec-02 sub-06 cat-01:
*    - addition max 64 bit; [sse2-avx512]; no epu
*  o sec-02 sub-06 cat-02:
*    - staturation addition; max 16 bit [sse2-avx512]; epu
*  o sec-02 sub-06 cat-03:
*    - subtraction; max 64bit [sse2-avx512]; no epu
*  o sec-02 sub-06 cat-04:
*    - staturation subtraction; max 16 bit; [sse2-avx512];u
*  o sec-02 sub-06 cat-05:
*    - Max
*  o sec-02 sub-06 cat-06:
*    - min
*  o sec-02 sub-06 cat-07:
*    - Horizontal max
*  o sec-02 sub-06 cat-08:
*    - Horizontal min
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-01:
+  - addition max 64 bit; [sse2-avx512]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define add_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm512_add_epi8((xVectI8), (yVectI8));)

#define add_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm512_add_epi16((xVectI16), (yVectI16));)

#define add_vectI32_retVectI32(xVectI32, yVectI32) \
   (_mm512_add_epi32((xVectI32), (yVectI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-02:
+  - staturation addition; max 16 bit [sse2 to avx512]; epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define addSat_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm512_adds_epu8((xVectU8), (yVectU8));)

#define addSat_vectU16_retVectU16(xVectU16, yVectU16)\
   (_mm512_adds_epu16((xVectU16), (yVectU16));)

#define addSat_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm512_adds_epi8((xVectI8), (yVectI8));)

#define addSat_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm512_adds_epi16((xVectI16), (yVectI16));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-03:
+  - subtraction; max 64bit [sse2-avx512]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sub_vectI8_retVectI8(xVectI8, yVectI8) \
    (_mm512_sub_epi8((xVectI8), (yVectI8));)

#define sub_vectI16_retVectI16(xVectI16, yVectI16) \
    (_mm512_sub_epi16((xVectI16), (yVectI16));)

#define sub_vectI32_retVectI32(xVectI32, yVectI32) \
    (_mm512_sub_epi32((xVectI32), (yVectI32));)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-04:
+  - staturation subtraction; max 16 bit; [sse2-avx512]); u
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define subSat_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm512_subs_epu8((xVectU8), (yVectU8));)

#define subSat_vectU16_retVectU16(xVectU16, yVectU16) \
   (_mm512_subs_epu16((xVectU16), (yVectU16));)

#define subSat_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm512_subs_epi8((xVectI8), (yVectI8));)

#define subSat_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm512_subs_epi16((xVectI16), (yVectI16));)

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-05:
+  - Max
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define max_vectI8_retVectI8(xVectI8,yVectI8)\
   (_mm512_max_epi8((xVectI8), (yVectI8));)

#define max_vectU8_retVectU8(xVectU8,yVectU8)\
   (_mm512_max_epu8((xVectU8), (yVectU8));)

#define max_vectI16_retVectI16(xVectI16,yVectI16)\
   (_mm512_max_epi16((xVectI16), (yVectI16));)

#define max_vectU16_retVectU16(xVectU16,yVectU16)\
   (_mm512_max_epu16((xVectU16), (yVectU16));)

#define max_vectI32_retVectI32(xVectI32,yVectI32)\
   (_mm512_max_epi32((xVectI32), (yVectI32));)

#define max_vectU32_retVectU32(xVectU32,yVectU32)\
   (_mm512_max_epu32((xVectU32), (yVectU32));)

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-06:
+  - Min
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define min_vectI8_retVectI8(xVectI8,yVectI8)\
   (_mm512_min_epi8((xVectI8), (yVectI8));)

#define min_vectU8_retVectU8(xVectU8,yVectU8)\
   (_mm512_min_epu8((xVectU8), (yVectU8));)

#define min_vectI16_retVectI16(xVectI16,yVectI16)\
   (_mm512_min_epi16((xVectI16), (yVectI16));)

#define min_vectU16_retVectU16(xVectU16,yVectU16)\
   (_mm512_min_epu16((xVectU16), (yVectU16));)

#define min_vectI32_retVectI32(xVectI32,yVectI32)\
   (_mm512_min_epi32((xVectI32), (yVectI32));)

#define min_vectU32_retVectU32(xVectU32,yVectU32)\
   (_mm512_min_epu32((xVectU32), (yVectU32));)

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-07:
+  - Horizontal max
+  - https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library 
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define hmax_vectI8_retI8(retI8, inVectI8) ({\
  vectI8 tmpVI8; /*To avoid overwriting input vector*/\
  \
  /*Find and move the max value to the lower part of the
  ` vector
  */\
  tmpVI8 =\
     max_vectI8_retVectI8(\
        (inVectI8), \
        srvect_vectI8_retVectI8((inVectI8), 32)\
  );\
  tmpVI8 =\
     max_vectI8_retVectI8(\
        (inVectI8), \
        srvect_vectI8_retVectI8((inVectI8), 16)\
  );\
  tmpVI8 =\
     max_vectI8_retVectI8(\
        (inVectI8), \
        srvect_vectI8_retVectI8((inVectI8), 8)\
  );\
  (tmpVI8) =\
     max_vectI8_retVectI8(\
        (tmpVI8), \
        srvect_vectI8_retVectI8((tmpVI8), 4)\
  );\
  (tmpVI8) =\
     max_vectI8_retVectI8(\
        (tmpVI8), \
        srvect_vectI8_retVectI8((tmpVI8), 2)\
  );\
  (tmpVI8) =\
     max_vectI8_retVectI8(\
        (tmpVI8), \
        srvect_vectI8_retVectI8((tmpVI8), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retI8 =\
     (int8_t) \
     extract_vectI16_I16(\
        cast_vectI8_to_VectI16((tmpVI8)),\
        0\
  );\
})

#define hmax_vectU8_retU8(retU8, inVectU8) ({\
  vectU8 tmpVU8; /*To avoid overwriting input vector*/\
  \
  /*Find and move the max value to the lower part of the
  ` vector
  */\
  tmpVU8 =\
     max_vectU8_retVectU8(\
        (inVectU8), \
        srvect_vectU8_retVectU8((inVectU8), 32)\
  );\
  tmpVU8 =\
     max_vectU8_retVectU8(\
        (inVectU8), \
        srvect_vectU8_retVectU8((inVectU8), 16)\
  );\
  tmpVU8 =\
     max_vectU8_retVectU8(\
        (inVectU8), \
        srvect_vectU8_retVectU8((inVectU8), 8)\
  );\
  (tmpVU8) =\
     max_vectU8_retVectU8(\
        (tmpVU8), \
        srvect_vectU8_retVectU8((tmpVU8), 4)\
  );\
  (tmpVU8) =\
     max_vectU8_retVectU8(\
        (tmpVU8), \
        srvect_vectU8_retVectU8((tmpVU8), 2)\
  );\
  (tmpVU8) =\
     max_vectU8_retVectU8(\
        (tmpVU8), \
        srvect_vectU8_retVectU8((tmpVU8), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retU8 =\
     (int8_t) \
     extract_vectI16_I16(\
        cast_vectU8_to_VectI16((tmpVU8)),\
        0\
  );\
})

#define hmax_vectI16_retI16(retI16, inVectI16) ({\
  vectI16 tmpVI16; /*To avoid overwriting input vector*/\
  \
  /*Find and move the max value to the lower part of the
  ` vector
  */\
  tmpVI16 =\
     max_vectI16_retVectI16(\
        (inVectI16), \
        srvect_vectI16_retVectI16((inVectI16), 16)\
  );\
  tmpVI16 =\
     max_vectI16_retVectI16(\
        (inVectI16), \
        srvect_vectI16_retVectI16((inVectI16), 8)\
  );\
  tmpVI16 =\
     max_vectI16_retVectI16(\
        (inVectI16), \
        srvect_vectI16_retVectI16((inVectI16), 4)\
  );\
  (tmpVI16) =\
     max_vectI16_retVectI16(\
        (tmpVI16), \
        srvect_vectI16_retVectI16((tmpVI16), 2)\
  );\
  (tmpVI16) =\
     max_vectI16_retVectI16(\
        (tmpVI16), \
        srvect_vectI16_retVectI16((tmpVI16), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retI16 = extract_vectI16_I16((tmpVI16), 0);\
}) /*hmaxVectI16*/

#define hmax_vectU16_retU16(retU16, inVectU16) ({\
  vectU16 tmpVU16; /*To avoid overwriting input vector*/\
  \
  /*Find and move the max value to the lower part of the
  ` vector
  */\
  tmpVU16 =\
     max_vectU16_retVectU16(\
        (inVectU16), \
        srvect_vectU16_retVectU16((inVectU16), 16)\
  );\
  tmpVU16 =\
     max_vectU16_retVectU16(\
        (inVectU16), \
        srvect_vectU16_retVectU16((inVectU16), 8)\
  );\
  tmpVU16 =\
     max_vectU16_retVectU16(\
        (inVectU16), \
        srvect_vectU16_retVectU16((inVectU16), 4)\
  );\
  (tmpVU16) =\
     max_vectU16_retVectU16(\
        (tmpVU16), \
        srvect_vectU16_retVectU16((tmpVU16), 2)\
  );\
  (tmpVU16) =\
     max_vectU16_retVectU16(\
        (tmpVU16), \
        srvect_vectU16_retVectU16((tmpVU16), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retU16 = extract_vectU16_U16((tmpVU16), 0);\
}) /*hmaxVectU16*/


/*Intel has a dedicated function for 32 bit, so I am
` assuming intel knows what they are doing here. It is
` based on their commands
*/

#define hmax_vectI32_retI32(retI32, inVectI32)\
   (_mm512_reduce_max_epi32((inVectI32)))

#define hmax_vectU32_retU32(retU32, inVectU32)\
   (_mm512_reduce_max_epu32((inVectU32)))

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-08:
+  - Horizontal min
+  - https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library 
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define hmin_vectI8_retI8(retI8, inVectI8) ({\
  vectI8 tmpVI8; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
  tmpVI8 =\
     min_vectI8_retVectI8(\
        (inVectI8), \
        srvect_vectI8_retVectI8((inVectI8), 32)\
  );\
  tmpVI8 =\
     min_vectI8_retVectI8(\
        (inVectI8), \
        srvect_vectI8_retVectI8((inVectI8), 16)\
  );\
  tmpVI8 =\
     min_vectI8_retVectI8(\
        (inVectI8), \
        srvect_vectI8_retVectI8((inVectI8), 8)\
  );\
  (tmpVI8) =\
     min_vectI8_retVectI8(\
        (tmpVI8), \
        srvect_vectI8_retVectI8((tmpVI8), 4)\
  );\
  (tmpVI8) =\
     min_vectI8_retVectI8(\
        (tmpVI8), \
        srvect_vectI8_retVectI8((tmpVI8), 2)\
  );\
  (tmpVI8) =\
     min_vectI8_retVectI8(\
        (tmpVI8), \
        srvect_vectI8_retVectI8((tmpVI8), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retI8 =\
     (int8_t) \
     extract_vectI16_I16(\
        cast_vectI8_to_VectI16((tmpVI8)),\
        0\
  );\
})

#define hmin_vectU8_retU8(retU8, inVectU8) ({\
  vectU8 tmpVU8; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
  tmpVU8 =\
     min_vectU8_retVectU8(\
        (inVectU8), \
        srvect_vectU8_retVectU8((inVectU8), 32)\
  );\
  tmpVU8 =\
     min_vectU8_retVectU8(\
        (inVectU8), \
        srvect_vectU8_retVectU8((inVectU8), 16)\
  );\
  tmpVU8 =\
     min_vectU8_retVectU8(\
        (inVectU8), \
        srvect_vectU8_retVectU8((inVectU8), 8)\
  );\
  (tmpVU8) =\
     min_vectU8_retVectU8(\
        (tmpVU8), \
        srvect_vectU8_retVectU8((tmpVU8), 4)\
  );\
  (tmpVU8) =\
     min_vectU8_retVectU8(\
        (tmpVU8), \
        srvect_vectU8_retVectU8((tmpVU8), 2)\
  );\
  (tmpVU8) =\
     min_vectU8_retVectU8(\
        (tmpVU8), \
        srvect_vectU8_retVectU8((tmpVU8), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retU8 =\
     (int8_t) \
     extract_vectI16_I16(\
        cast_vectU8_to_VectI16((tmpVU8)),\
        0\
  );\
})

#define hmin_vectI16_retI16(retI16, inVectI16) ({\
  vectI16 tmpVI16; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
  tmpVI16 =\
     min_vectI16_retVectI16(\
        (inVectI16), \
        srvect_vectI16_retVectI16((inVectI16), 16)\
  );\
  tmpVI16 =\
     min_vectI16_retVectI16(\
        (inVectI16), \
        srvect_vectI16_retVectI16((inVectI16), 8)\
  );\
  tmpVI16 =\
     min_vectI16_retVectI16(\
        (inVectI16), \
        srvect_vectI16_retVectI16((inVectI16), 4)\
  );\
  (tmpVI16) =\
     min_vectI16_retVectI16(\
        (tmpVI16), \
        srvect_vectI16_retVectI16((tmpVI16), 2)\
  );\
  (tmpVI16) =\
     min_vectI16_retVectI16(\
        (tmpVI16), \
        srvect_vectI16_retVectI16((tmpVI16), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retI16 = extract_vectI16_I16((tmpVI16)), 0);\
}) /*hminVectI16*/

#define hmin_vectU32_retU32(retU32, inVectU32) ({\
  vectU32 tmpVU32; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
  tmpVU32 =\
     min_vectU32_retVectU32(\
        (inVectU32), \
        srvect_vectU32_retVectU32((inVectU32), 16)\
  );\
  tmpVU32 =\
     min_vectU32_retVectU32(\
        (inVectU32), \
        srvect_vectU32_retVectU32((inVectU32), 8)\
  );\
  tmpVU32 =\
     min_vectU32_retVectU32(\
        (inVectU32), \
        srvect_vectU32_retVectU32((inVectU32), 4)\
  );\
  (tmpVU32) =\
     min_vectU32_retVectU32(\
        (tmpVU32), \
        srvect_vectU32_retVectU32((tmpVU32), 2)\
  );\
  (tmpVU32) =\
     min_vectU32_retVectU32(\
        (tmpVU32), \
        srvect_vectU32_retVectU32((tmpVU32), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retU32 = extract_vectU32_U32((tmpVU32)), 0);\
}) /*hminVectU32*/

/*Intel has a dedicated function for 32 bit, so I am
` assuming intel knows what they are doing here.
*/

#define hmin_vectI32_retI32(retI32, inVectI32)\
   (_mm512_reduce_min_epi32((inVectI32)))

#define hmin_vectU32_retU32(retU32, inVectU32)\
   (_mm512_reduce_min_epu32((inVectU32)))


/*********************************************************\
* Sec-02 Sub-07:
*  - Conversions [Here for NEON support]
*  o sec02 sub-07 cat-01:
*    - Convert signed vectors to other signed vectors
*  o sec02 sub-07 cat-02:
*    - Convert unsigned vectors to other unsigned vectors
*  o sec02 sub-07 cat-03:
*    - Convert unsigned vectors to signed vectors
*  o sec02 sub-07 cat-04:
*    - Convert signed vectors to unsigned vectors
*  o sec02 sub-07 cat-05:
*    - Convert masks to vectors
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-01:
+  - Convert signed vectors to other signed vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectI8_to_VectI16(inVectI8) \
  ((vectI16) (inVectI8);)

#define cast_vectI8_to_VectI32(inVectI8) \
  ((vectI32) (inVectI8);)

#define cast_vectI16_to_VectI8(inVectI16) \
  ((vectI8) (inVectI16);)

#define cast_vectI16_to_VectI32(inVectI16) \
  ((vectI32) (inVectI16);)

#define cast_vectI32_to_VectI8(inVectI32) \
  ((vectI8) (inVectI32);)

#define cast_vectI32_to_VectI16(inVectI32) \
  ((vectI16) (inVectI32);)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-02:
+  - Convert unsigned vectors to unsigned vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectU8_to_VectU16(inVectU8) \
  ((vectU16) (inVectU8);)

#define cast_vectU8_to_VectU32(inVectU8) \
  ((vectU32) (inVectU8);)

#define cast_vectU16_to_VectU8(inVectU16) \
  ((vectU8) (inVectU16);)

#define cast_vectU16_to_VectU32(inVectU16) \
  ((vectU32) (inVectU16);)

#define cast_vectU32_to_VectU8(inVectU32) \
  ((vectU8) (inVectU32);)

#define cast_vectU32_to_VectU16(inVectU32) \
  ((vectU16) (inVectU32);)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-03:
+  - Convert unsigned vectors to signed vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectU8_to_VectI8(inVectU8) \
  ((vectI8) (inVectU8);)

#define cast_vectU8_to_VectI16(inVectU8) \
  ((vectI16) (inVectU8);)

#define cast_vectU8_to_VectI32(inVectU8) \
  ((vectI32) (inVectU8);)

#define cast_vectU16_to_VectI8(inVectU16) \
  ((vectI8) (inVectU16);)

#define cast_vectU16_to_VectI16(inVectU16) \
  ((vectI16) (inVectU16);)

#define cast_vectU16_to_VectI32(inVectU16) \
  ((vectI32) (inVectU16);)

#define cast_vectU32_to_VectI8(inVectU32) \
  ((vectI8) (inVectU32);)

#define cast_vectU32_to_VectI16(inVectU32) \
  ((vectI16) (inVectU32);)

#define cast_vectU32_to_VectI32(inVectU32) \
  ((vectI32) (inVectU32);)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-04:
+  - Convert signed vectors to unsigned vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectI8_to_VectU8(inVectI8) \
  ((vectU8) (inVectI8);)

#define cast_vectI8_to_VectU16(inVectI8) \
  ((vectU16) (inVectI8);)

#define cast_vectI8_to_VectU32(inVectI8) \
  ((vectU32) (inVectI8);)

#define cast_vectI16_to_VectU8(inVectI16) \
  ((vectU8) (inVectI16);)

#define cast_vectI16_to_VectU16(inVectI16) \
  ((vectU16) (inVectI16);)

#define cast_vectI16_to_VectU32(inVectI16) \
  ((vectU32) (inVectI16);)

#define cast_vectI32_to_VectU8(inVectI32) \
  ((vectU8) (inVectI32);)

#define cast_vectI32_to_VectU16(inVectI32) \
  ((vectU16) (inVectI32);)

#define cast_vectI32_to_VectU32(inVectI32) \
  ((vectU32) (inVectI32);)

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-05:
+  - Convert masks to vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define mask8_retVectI8(inMask8)\
   ((vectI8) _mm512_movm_epi8((inMask8));)

#define mask16_retVectI16(inMask16)\
   ((vectI16) _mm512_movm_epi16((inMask16));)

#define mask32_retVectI32(inMask32)\
   ((vectI32) _mm512_movm_epi32((inMask32));)

#define mask8_retVectU8(inMask8)\
   ((vectU8) _mm512_movm_epi8((inMask8));)

#define mask16_retVectU16(inMask16)\
   ((vectU16) _mm512_movm_epi16((inMask16));)

#define mask32_retVectU32(inMask32)\
   ((vectU32) _mm512_movm_epi32((inMask32));)

#endif
