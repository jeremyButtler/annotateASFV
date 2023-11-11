/*#########################################################
# Name: vecWrapAvx2Fun.h
# Use:
#  - Holds macro wrappers and variables for AVX2 intrinsics
# Note:
#  - Do not use the min/max for unsigned 32 ints if you
#    plan on the user using SSE2. At best it is as efficent
#    as doing it without vectors. With loading costs, it
#    is worse.
#########################################################*/

#ifndef VECTWRAPAVX2_H
#define VECTWRAPAVX2_H 

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' vectWrapAvx2Fun.h TOC: Table Of Contents
'  - Holds macro wrapper for AVX2 intrinsics
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

#define defVectBits 256
#define defVectBytes 32
#define defNum8BitElms 32
#define defNum16BitElms 16
#define defNum32BitElms 8
#define defMvMask8 0xffffffff
#define defMvMask16 0xffffffff
#define defMvMask32 0xffffffff

typedef __m256i vectI8;  /*vector of bytes (8bits)*/
typedef __m256i vectI16; /*vector of shorts (16 bits)*/
typedef __m256i vectI32; /*vector of ints (32 bits)*/

typedef __m256i vectU8;  /*Prefer unsigned bytes*/
typedef __m256i vectU16; /*Prefer unsigned shorts*/
typedef __m256i vectU32; /*Prefer unsinged ints*/

typedef __m256i mask8;  /*mask of 8 bit values*/
typedef __m256i mask16; /*Mask of 16 bit values*/
typedef __m256i mask32; /*mask of 32 bit values*/

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
  (_mm256_load_si256((__m256i *) (arrayI8)))

#define load_I16Ary_retVectI16(arrayI16) \
  (_mm256_load_si256((__m256i *) (arrayI16)))

#define load_I32Ary_retVectI32(arrayI32) \
  (_mm256_load_si256((__m256i *) (arrayI32)))


#define load_U8Ary_retVectU8(arrayU8) \
  (_mm256_load_si256((__m256i *) (arrayU8)))

#define load_U16Ary_retVectU16(arrayU16) \
  (_mm256_load_si256((__m256i *) (arrayU16)))

#define load_U32Ary_retVectU32(arrayU32) \
  (_mm256_load_si256((__m256i *) (arrayU32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-02:
+  - Unaligned loading
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#define loadu_I8Ary_retVectI8(arrayI8) \
  (_mm256_loadu_si256((__m256i *) (arrayI8)))

#define loadu_I16Ary_retVectI16(arrayI16) \
  (_mm256_loadu_si256((__m256i *) (arrayI16)))

#define loadu_I32Ary_retVectI32(arrayI32) \
  (_mm256_loadu_si256((__m256i *) (arrayI32)))


#define loadu_U8Ary_retVectU8(arrayU8) \
  (_mm256_loadu_si256((__m256i *) (arrayU8)))

#define loadu_U16Ary_retVectU16(arrayU16) \
  (_mm256_loadu_si256((__m256i *) (arrayU16)))

#define loadu_U32Ary_retVectU32(arrayU32) \
  (_mm256_loadu_si256((__m256i *) (arrayU32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-03:
+  - make zero vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define zero_retVectI8()  (_mm256_setzero_si256())
#define zero_retVectI16() (_mm256_setzero_si256())
#define zero_retVectI32() (_mm256_setzero_si256())

#define zero_retVectU8()  (_mm256_setzero_si256())
#define zero_retVectU16() (_mm256_setzero_si256())
#define zero_retVectU32() (_mm256_setzero_si256())

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-04:
+  - Make vectors of one element
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define set1_I8_retVectI8(valC) (_mm256_set1_epi8((valC)))

#define set1_I16_retVectI16(valI16) \
   (_mm256_set1_epi16((valI16)))

#define set1_I32_retVectI32(valI32) \
   (_mm256_set1_epi32((valI32)))


#define set1_U8_retVectU8(valUC) \
   (_mm256_set1_epi8((valUC)))

#define set1_U16_retVectU16(valC) \
   (_mm256_set1_epi16((valU16)))

#define set1_U32_retVectU32(valU32) \
   (_mm256_set1_epi32((valU32)))


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-01 Cat-05:
+  - Insert an element into an vector
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define insert_I16_retVectI16(insVectI16, insI16, posI)\
   (_mm256_insert_epi16((insVectI16), (insI16), (posI)))

#define insert_I32_retVectI32(insVectI32,insI32,posI)\
   (_mm256_insert_epi32((insVectI32), (insI32), (posI)))

#define insert_U16_retVectU16(insVectU16, insU16, posU)\
   (_mm256_insert_epi16((insVectU16), (insU16), (posU)))

#define insert_U32_retVectU32(insVectU32,insU32,posU)\
   (_mm256_insert_epi32((insVectU32), (insU32), (posU)))

/*********************************************************\
* Sec-02 Sub-02:
*  - Output functions
*  o sec-02 sub-02 cat-01:
*    - Store vector ouput into an array
*  o sec-02 sub-02 cat-02:
*    - Store vector ouput into an array (unaligned)
*  o sec-02 sub-02 cat-03:
*    - Store masks into longs
*  o sec-02 sub-02 cat-04:
*    - Extract an element
\*********************************************************/

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-01:
+  - Store vector ouput into an array (aligned)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define store_vectI8_retAryI8(retAry, inVectI8) \
  (_mm256_store_si256((__m256i *) (retAry), (inVectI8)))

#define store_vectI16_retAryI16(retAry, inVectI16) \
  (_mm256_store_si256((__m256i *) (retAry), (inVectI16)))

#define store_vectI32_retAryI32(retAry, inVectI32) \
  (_mm256_store_si256((__m256i *) (retAry), (inVectI32)))

#define store_vectU8_retAryU8(retAry, inVectU8) \
  (_mm256_store_si256((__m256i *) (retAry), (inVectU8)))

#define store_vectU16_retAryU16(retAry, inVectU16) \
  (_mm256_store_si256((__m256i *) (retAry), (inVectU16)))

#define store_vectU32_retAryU32(retAry, inVectU32) \
  (_mm256_store_si256((__m256i *) (retAry), (inVectU32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-02:
+  - Store vector ouput into an array (unaligned)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define storeu_vectI8_retAryI8(retAry, inVectI8) \
  (_mm256_storeu_si256((__m256i *) (retAry), (inVectI8)))

#define storeu_vectI16_retAryI16(retAry, inVectI16) \
  (_mm256_storeu_si256((__m256i *) (retAry), (inVectI16)))

#define storeu_vectI32_retAryI32(retAry, inVectI32) \
  (_mm256_storeu_si256((__m256i *) (retAry), (inVectI32)))

#define storeu_vectU8_retAryU8(retAry, inVectU8) \
  (_mm256_storeu_si256((__m256i *) (retAry), (inVectU8)))

#define storeu_vectU16_retAryU16(retAry, inVectU16) \
  (_mm256_storeu_si256((__m256i *) (retAry), (inVectU16)))

#define storeu_vectU32_retAryU32(retAry, inVectU32) \
  (_mm256_storeu_si256((__m256i *) (retAry), (inVectU32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-03:
+  - Store masks into longs
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*storing in retUL because this is how AVX512 does it*/

#define store_mask8_retUL(retUL, mask8) \
   ((retUL) = _mm256_movemask_epi8((mask8)))

#define store_mask16_retUL(retUL, mask16) \
   ((retUL) = _mm256_movemask_epi8((mask16)))

#define store_mask32_retUL(retUL, mask32) \
   ((retUL) = _mm256_movemask_epi8((mask32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-02 Cat-04:
+  - Extract an element
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define extract_vectI8_I8(inVectI8, posI)\
   (_mm256_extract_epi8((inVectI8), (posI)))

#define extract_vectU8_U8(inVectU8, posI)\
   (_mm256_extract_epi8((inVectU8), (posI)))

#define extract_vectI16_I16(inVectI16, posI)\
   (_mm256_extract_epi16((inVectI16), (posI)))

#define extract_vectU16_U16(inVectU16, posI)\
   (_mm256_extract_epi16((inVectU16), (posI)))

#define extract_vectI32_I32(inVectI32, posI)\
   (_mm256_extract_epi32((inVectI32), (posI)))

#define extract_vectU32_U32(inVectU32, posI)\
   (_mm256_extract_epi32((inVectU32), (posI)))

/*********************************************************\
* Sec-02 Sub-03:
*  - Comparison functions (if)
*  o sec-02 sub-03 cat-01:
*    - Comparisons (return masks)
*  o sec-02 sub-03 cat-02:
*    - Unsinged equal comparisions (returns masks)
*  o sec-02 sub-03 cat-03:
*    - Unsinged less than comparisions (less efficentmask)
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

/*8-bit comparisions returns a mask8*/

#define cmpeq_vectI8_retMask8(xVectI8, yVectI8) \
   ((mask8) _mm256_cmpeq_epi8((xVectI8), (yVectI8)))

#define cmpgt_vectI8_retMask8(xVectI8, yVectI8) \
   ((mask8) _mm256_cmpgt_epi8((xVectI8), (yVectI8)))

#define cmplt_vectI8_retMask8(xVectI8, yVectI8) \
   ((mask8) _mm256_cmplt_epi8((xVectI8), (yVectI8)))

/*16-bit comparisions returns a mask16*/

#define cmpeq_vectI16_retMask16(xVectI16, yVectI16) \
  ((mask16) _mm256_cmpeq_epi16((xVectI16), (yVectI16)))

#define cmpgt_vectI16_retMask16(xVectI16, yVectI16) \
  ((mask16) _mm256_cmpgt_epi16((xVectI16), (yVectI16)))

#define cmplt_vectI16_retMask16(xVectI16, yVectI16) \
  ((mask16) _mm256_cmplt_epi16((xVectI16), (yVectI16)))

/*32-bit comparisions returns a mask32*/

#define cmpeq_vectI32_retMask32(xVectI32, yVectI32) \
  ((mask32) _mm256_cmpeq_epi32((xVectI32), (yVectI32)))

#define cmpgt_vectI32_retMask32(xVectI32, yVectI32) \
  ((mask32) _mm256_cmpgt_epi32((xVectI32), (yVectI32)))

#define cmplt_vectI32_retMask32(xVectI32, yVectI32) \
  ((mask32) _mm256_cmplt_epi32((xVectI32), (yVectI32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-02:
+  - Unsinged equal comparisions (returns masks)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*These are the one to one equivlent comparisons*/

#define cmpeq_vectU8_retMask8(xVectI8, yVectI8) \
   ((mask8) _mm256_cmpeq_epi8((xVectI8), (yVectI8)))

#define cmpeq_vectU16_retMask16(xVectI16, yVectI16) \
  ((mask16) _mm256_cmpeq_epi16((xVectI16), (yVectI16)))

#define cmpeq_vectU32_retMask32(xVectI32, yVectI32) \
  ((mask32) _mm256_cmpeq_epi32((xVectI32), (yVectI32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-03:
+  - Unsinged less than comparisions (less efficent; mask)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*< comparisons require 5 op & 3 extra vectors.*/

#define cmplt_vectU8_retMask8(xVectU8, yVectU8)({\
   vectU8 ltCmpV = _mm256_cmplt_epi8(xVectU8, yVectU8);\
     /*-x < -y when unsigned x > unsinged y
     `  x <  y when unsinged x < unsigned y
     `  x > -y when unsigned x < unsinged y
     ` -x <  y when unsigned x > unsigned y
     */\
   vectU8 zeroV = zero_retVectU8();\
   \
   vectU8 negCmpV = _mm256_cmplt_epi8(xVectU8, zeroV);\
   zeroV = _mm256_cmplt_epi8(yVectU8, zeroV);\
   negCmpV = _mm256_or_si256(negCmpV, zeroV);\
      /*Captures -1 if x or y is negative; 0 otherwise*/\
   _mm256_xor_si256(ltCmpV, negCmpV);\
      /*Results go from -1 to 0 or 0 to -1 when x or y is
      ` negative, but stay the same if x and y are positive
      */\
})

#define cmplt_vectU16_retMask16(xVectU16, yVectU16)({\
   vectU16 ltCmpV = _mm256_cmplt_epi16(xVectU16,yVectU16);\
     /*-x < -y when unsigned x > unsinged y
     `  x <  y when unsinged x < unsigned y
     `  x > -y when unsigned x < unsinged y
     ` -x <  y when unsigned x > unsigned y
     */\
   vectU16 zeroV = zero_retVectU16();\
   \
   vectU16 negCmpV = _mm256_cmplt_epi16(xVectU16, zeroV);\
   zeroV = _mm256_cmplt_epi16(yVectU16, zeroV);\
   negCmpV = _mm256_or_si256(negCmpV, zeroV);\
      /*Captures -1 if x or y is negative; 0 otherwise*/\
   _mm256_xor_si256(ltCmpV, negCmpV);\
      /*Results go from -1 to 0 or 0 to -1 when x or y is
      ` negative, but stay the same if x and y are positive
      */\
})

/*See cmplt_vectU8_retMask8 for logic*/
#define cmplt_vectU32_retMask32(xVectU32, yVectU32)({\
   vectU32 ltCmpV = _mm256_cmplt_epi32(xVectU32,yVectU32);\
     /*-x < -y when unsigned x > unsinged y
     `  x <  y when unsinged x < unsigned y
     `  x > -y when unsigned x < unsinged y
     ` -x <  y when unsigned x > unsigned y
     */\
   vectU32 zeroV = zero_retVectU32();\
   \
   vectU32 negCmpV = _mm256_cmplt_epi32(xVectU32, zeroV);\
   zeroV = _mm256_cmplt_epi32(yVectU32, zeroV);\
   negCmpV = _mm256_or_si256(negCmpV, zeroV);\
      /*Captures -1 if x or y is negative; 0 otherwise*/\
   _mm256_xor_si256(ltCmpV, negCmpV);\
      /*Results go from -1 to 0 or 0 to -1 when x or y is
      ` negative, but stay the same if x and y are positive
      */\
})

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-04:
+  - Unsinged greater than comparisions (less efficent)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/*> comparisons require 5 op & 3 extra vectors.*/

#define cmpgt_vectU8_retMask8(xVectU8, yVectU8)({\
   vectU8 gtCmpV = _mm256_cmpgt_epi8(xVectU8, yVectU8);\
     /*-x > -y when unsigned x < unsinged y
     `  x >  y when unsinged x > unsigned y
     `  x <  y when unsinged x < unsigned y
     `  x > -y when unsigned x < unsinged y
     ` -x <  y when unsigned x > unsigned y
     */\
   vectU8 zeroV = zero_retVectU8();\
   \
   vectU8 negCmpV = _mm256_cmplt_epi8(xVectU8, zeroV);\
   zeroV = _mm256_cmplt_epi8(yVectU8, zeroV);\
   negCmpV = _mm256_or_si256(negCmpV, zeroV);\
      /*Captures -1 if x or y is negative; 0 otherwise*/\
   _mm256_xor_si256(gtCmpV, negCmpV);\
      /*Results go from -1 to 0 or 0 to -1 when x or y is
      ` negative, but stay the same if x and y are positive
      */\
})

#define cmpgt_vectU16_retMask16(xVectU16, yVectU16)({\
   vectU16 gtCmpV = _mm256_cmpgt_epi16(xVectU16,yVectU16);\
     /*-x < -y when unsigned x > unsinged y
     `  x >  y when unsinged x > unsigned y
     `  x <  y when unsinged x < unsigned y
     `  x > -y when unsigned x < unsinged y
     ` -x <  y when unsigned x > unsigned y
     */\
   vectU16 zeroV = zero_retVectU16();\
   \
   vectU16 negCmpV = _mm256_cmplt_epi16(xVectU16, zeroV);\
   zeroV = _mm256_cmplt_epi16(yVectU16, zeroV);\
   negCmpV = _mm256_or_si256(negCmpV, zeroV);\
      /*Captures -1 if x or y is negative; 0 otherwise*/\
   _mm256_xor_si256(gtCmpV, negCmpV);\
      /*Results go from -1 to 0 or 0 to -1 when x or y is
      ` negative, but stay the same if x and y are positive
      */\
})

/*See cmplt_vectU8_retMask8 for logic*/
#define cmpgt_vectU32_retMask32(xVectU32, yVectU32)({\
   vectU32 gtCmpV = _mm256_cmpgt_epi32(xVectU32,yVectU32);\
     /*-x < -y when unsigned x > unsinged y
     `  x >  y when unsinged x > unsigned y
     `  x <  y when unsinged x < unsigned y
     `  x > -y when unsigned x < unsinged y
     ` -x <  y when unsigned x > unsigned y
     */\
   vectU32 zeroV = zero_retVectU32();\
   \
   vectU32 negCmpV = _mm256_cmplt_epi32(xVectU32, zeroV);\
   zeroV = _mm256_cmplt_epi32(yVectU32, zeroV);\
   negCmpV = _mm256_or_si256(negCmpV, zeroV);\
      /*Captures -1 if x or y is negative; 0 otherwise*/\
   _mm256_xor_si256(gtCmpV, negCmpV);\
      /*Results go from -1 to 0 or 0 to -1 when x or y is
      ` negative, but stay the same if x and y are positive
      */\
})

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-03 Cat-05:
+  - Fix differences in population counts (total 1's)
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/* This is here to support SSE2, AVX2, and NEON. 
`  The returned mask in AVX128 has only one bit per
`  data type in comparison
`  SSE (all) and AVX2 returns a 128 bit vector that is
`   converted to characters
*/
#define fix_mask8_popcount(inUL)  ((inUL))
#define fix_mask16_popcount(inUL) ((inUL) >> 1)
#define fix_mask32_popcount(inUL) ((inUL) >> 2)

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
   (_mm256_andnot_si256((notVectI8), (inVectI8)))

#define andNot_vectU8_retVectU8(notVectU8, inVectU8) \
   (_mm256_andnot_si256((notVectU8), (inVectU8)))

#define andNot_vectI16_retVectI16(notVectI16, inVectI16) \
   (_mm256_andnot_si256((notVectI16), (inVectI16)))

#define andNot_vectU16_retVectU16(notVectU16, inVectU16) \
   (_mm256_andnot_si256((notVectU16), (inVectU16)))

#define andNot_vectI32_retVectI32(notVectI32, inVectI32) \
   (_mm256_andnot_si256((notVectI32), (inVectI32)))

#define andNot_vectU32_retVectU32(notVectU32, inVectU32) \
   (_mm256_andnot_si256((notVectU32), (inVectU32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-02:
+  - and functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define and_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm256_and_si256((xVectI8), (yVectI8)))

#define and_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm256_and_si256((xVectU8), (yVectU8)))

#define and_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm256_and_si256((xVectI16), (yVectI16)))

#define and_vectU16_retVectU16(xVectU16, yVectU16) \
   (_mm256_and_si256((xVectU16), (yVectU16)))

#define and_vectI32_retVectI32(xVectI32, yVectI32) \
   (_mm256_and_si256((xVectI32), (yVectI32)))

#define and_vectU32_retVectU32(xVectU32, yVectU32) \
   (_mm256_and_si256((xVectU32), (yVectU32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-03:
+  - or functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define or_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm256_or_si256((xVectI8), (yVectI8)))

#define or_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm256_or_si256((xVectU8), (yVectU8)))

#define or_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm256_or_si256((xVectI16), (yVectI16)))

#define or_vectU16_retVectU16(xVectU16, yVectU16) \
   (_mm256_or_si256((xVectU16), (yVectU16)))

#define or_vectI32_retVectI32(xVectI32, yVectI32) \
   (_mm256_or_si256((xVectI32), (yVectI32)))

#define or_vectU32_retVectU32(xVectU32, yVectU32) \
   (_mm256_or_si256((xVectU32), (yVectU32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-04 Cat-04:
+  - xor functions
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define xor_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm256_xor_si256((xVectI8), (yVectI8)))

#define xor_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm256_xor_si256((xVectU8), (yVectU8)))

#define xor_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm256_xor_si256((xVectI16), (yVectI16)))

#define xor_vectU16_retVectU16(xVectU16, yVectU16) \
   (_mm256_xor_si256((xVectU16), (yVectU16)))

#define xor_vectI32_retVectI32(xVectI32, yVectI32) \
   (_mm256_xor_si256((xVectI32), (yVectI32)))

#define xor_vectU32_retVectU32(xVectU32, yVectU32) \
   (_mm256_xor_si256((xVectU32), (yVectU32)))

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
    (_mm256_srli_epi8((inVectI8), (numBitsI)))

#define sr_vectU8_retVectU8(inVectU8, numBitsI)\
    (_mm256_srli_epi8((inVectU8), (numBitsI)))

#define sr_vectI16_retVectI16(inVectI16, numBitsI)\
    (_mm256_srli_epi16((inVectI16), (numBitsI)))

#define sr_vectU16_retVectU16(inVectU16, numBitsI)\
    (_mm256_srli_epi16((inVectU16), (numBitsI)))

#define sr_vectI32_retVectI32(inVectI32, numBitsI)\
    (_mm256_srli_epi32((inVectI32), (numBitsI)))

#define sr_vectU32_retVectU32(inVectU32, numBitsI)\
    (_mm256_srli_epi32((inVectU32), (numBitsI)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-02:
+  - Shift a vector right by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define srvect_vectI8_retVectI8(inVectI8, numBytesI)\
    (_mm256_srli_si256((inVectI8), (numBytesI)))

#define srvect_vectU8_retVectU8(inVectU8, numBytesI)\
    (_mm256_srli_si256((inVectU8), (numBytesI)))

#define srvect_vectI16_retVectI16(inVectI16, numShortsI)\
    (_mm256_srli_si256((inVectI16), (numShortsI) << 1))

#define srvect_vectU16_retVectU16(inVectU16, numShortsI)\
    (_mm256_srli_si256((inVectU16), (numShortsI) << 1))

#define srvect_vectI32_retVectI32(inVectI32, numIntsI)\
    (_mm256_srli_si256((inVectI32), (numIntsI) << 2))

#define srvect_vectU32_retVectU32(inVectU32, numIntsI)\
    (_mm256_srli_si256((inVectU32), (numIntsI) << 2))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-03:
+  - Shift numbers in vectors left by x bits
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sl_vectI8_retVectI8(inVectI8, numBitsI)\
    (_mm256_slli_epi8((inVectI8), (numBitsI)))

#define sl_vectU8_retVectU8(inVectU8, numBitsI)\
    (_mm256_slli_epi8((inVectU8), (numBitsI)))

#define sl_vectI16_retVectI16(inVectI16, numBitsI)\
    (_mm256_slli_epi16((inVectI16), (numBitsI)))

#define sl_vectU16_retVectU16(inVectU16, numBitsI)\
    (_mm256_slli_epi16((inVectU16), (numBitsI)))

#define sl_vectI32_retVectI32(inVectI32, numBitsI)\
    (_mm256_slli_epi32((inVectI32), (numBitsI)))

#define sl_vectU32_retVectU32(inVectU32, numBitsI)\
    (_mm256_slli_epi32((inVectU32), (numBitsI)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-05 Cat-04:
+  - Shift a vector left by x bytes
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define slvect_vectI8_retVectI8(inVectI8, numBytesI)\
    (_mm256_slli_si256((inVectI8), (numBytesI)))

#define slvect_vectU8_retVectU8(inVectU8, numBytesI)\
    (_mm256_slli_si256((inVectU8), (numBytesI)))

#define slvect_vectI16_retVectI16(inVectI16, numShortsI)\
    (_mm256_slli_si256((inVectI16), (numShortsI) << 1))

#define slvect_vectU16_retVectU16(inVectU16, numShortsI)\
    (_mm256_slli_si256((inVectU16), (numShortsI) << 1))

#define slvect_vectI32_retVectI32(inVectI32, numIntsI)\
    (_mm256_slli_si256((inVectI32), (numIntsI) << 2))

#define slvect_vectU32_retVectU32(inVectU32, numIntsI)\
    (_mm256_slli_si256((inVectU32), (numIntsI) << 2))

/*********************************************************\
* Sec-02 Sub-06:
*  - Math functions
*  o sec-02 sub-06 cat-01:
*    - addition max 64 bit [sse2-avx128]; no epu
*  o sec-02 sub-06 cat-02:
*    - staturation addition max 16 bit [sse2-avx128]; epu
*  o sec-02 sub-06 cat-03:
*    - subtraction max 64bit [sse2-avx128]; no epu
*  o sec-02 sub-06 cat-04:
*    - staturation subtraction max 16 bit; [sse2-avx128];u
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
+  - addition max 64 bit [sse2-avx128]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define add_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm256_add_epi8((xVectI8), (yVectI8)))

#define add_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm256_add_epi16((xVectI16), (yVectI16)))

#define add_vectI32_retVectI32(xVectI32, yVectI32) \
   (_mm256_add_epi32((xVectI32), (yVectI32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-02:
+  - staturation addition max 16 bit [sse2 to avx128]; epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define addSat_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm256_adds_epu8((xVectU8), (yVectU8)))

#define addSat_vectU16_retVectU16(xVectU16, yVectU16)\
   (_mm256_adds_epu16((xVectU16), (yVectU16)))

#define addSat_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm256_adds_epi8((xVectI8), (yVectI8)))

#define addSat_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm256_adds_epi16((xVectI16), (yVectI16)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-03:
+  - subtraction max 64bit [sse2-avx128]; no epu
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define sub_vectI8_retVectI8(xVectI8, yVectI8) \
    (_mm256_sub_epi8((xVectI8), (yVectI8)))

#define sub_vectI16_retVectI16(xVectI16, yVectI16) \
    (_mm256_sub_epi16((xVectI16), (yVectI16)))

#define sub_vectI32_retVectI32(xVectI32, yVectI32) \
    (_mm256_sub_epi32((xVectI32), (yVectI32)))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-04:
+  - staturation subtraction max 16 bit; [sse2-avx128]); u
\++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define subSat_vectU8_retVectU8(xVectU8, yVectU8) \
   (_mm256_subs_epu8((xVectU8), (yVectU8)))

#define subSat_vectU16_retVectU16(xVectU16, yVectU16) \
   (_mm256_subs_epu16((xVectU16), (yVectU16)))

#define subSat_vectI8_retVectI8(xVectI8, yVectI8) \
   (_mm256_subs_epi8((xVectI8), (yVectI8)))

#define subSat_vectI16_retVectI16(xVectI16, yVectI16) \
   (_mm256_subs_epi16((xVectI16), (yVectI16)))

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-05:
+  - Max
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define max_vectI8_retVectI8(xVectI8,yVectI8)\
   (_mm256_max_epi8((xVectI8), (yVectI8)))

#define max_vectU8_retVectU8(xVectU8,yVectU8)\
   (_mm256_max_epu8((xVectU8), (yVectU8)))

#define max_vectI16_retVectI16(xVectI16,yVectI16)\
   (_mm256_max_epi16((xVectI16), (yVectI16)))

#define max_vectU16_retVectU16(xVectU16,yVectU16)\
   (_mm256_max_epu16((xVectU16), (yVectU16)))

#define max_vectI32_retVectI32(xVectI32,yVectI32)\
   (_mm256_max_epi32((xVectI32), (yVectI32)))

#define max_vectU32_retVectU32(xVectU32,yVectU32)\
   (_mm256_max_epu32((xVectU32), (yVectU32)))

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-06:
+  - Min
\+++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define min_vectI8_retVectI8(xVectI8,yVectI8)\
  (_mm256_min_epi8((xVectI8), (yVectI8)))

#define min_vectU8_retVectU8(xVectU8,yVectU8)\
   (_mm256_min_epu8((xVectU8), (yVectU8)))

#define min_vectI16_retVectI16(xVectI16,yVectI16)\
   (_mm256_min_epi16((xVectI16), (yVectI16)))

#define min_vectU16_retVectU16(xVectU16,yVectU16)\
   (_mm256_min_epu16((xVectU16), (yVectU16)))

#define min_vectI32_retVectI32(xVectI32,yVectI32)\
   (_mm256_min_epi32((xVectI32), (yVectI32)))

#define min_vectU32_retVectU32(xVectU32,yVectU32)\
   (_mm256_min_epu32((xVectU32), (yVectU32)))

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
}) /*hmaxVectI8*/


#define hmax_vectU8_retU8(retU8, inVectU8) ({\
  vectU8 tmpVU8; /*To avoid overwriting input vector*/\
  \
  /*Find and move the max value to the lower part of the
  ` vector
  */\
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
     (uint8_t) \
     extract_vectI16_I16(\
        cast_vectU8_to_VectI16((tmpVU8)),\
        0\
  );\
}) /*hmaxVectU8*/

#define hmax_vectI16_retI16(retI16, inVectI16) ({\
  vectI16 tmpVI16; /*To avoid overwriting input vector*/\
  \
  /*Find and move the max value to the lower part of the
  ` vector
  */\
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

#define hmax_vectI32_retI32(retI32, inVectI32) ({\
  vectI32 tmpVI32; /*To avoid overwriting input vector*/\
  \
  /*Find and move the max value to the lower part of the
  ` vector
  */\
  tmpVI32 =\
     max_vectI32_retVectI32(\
        (inVectI32), \
        srvect_vectI32_retVectI32((inVectI32), 4)\
  );\
  tmpVI32 =\
     max_vectI32_retVectI32(\
        (inVectI32), \
        srvect_vectI32_retVectI32((inVectI32), 2)\
  );\
  (tmpVI32) =\
     max_vectI32_retVectI32(\
        (tmpVI32), \
        srvect_vectI32_retVectI32((tmpVI32), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retI32 = extract_vectI32_I32((tmpVI32), 0);\
   /*extract 32I is a double extract for SSE2 (not SSE4)*/\
}) /*hmaxVectI32*/

#define hmax_vectU32_retU32(retU32, inVectU32) ({\
  vectU32 tmpVU32; /*To avoid overwriting input vector*/\
  \
  /*Find and move the max value to the lower part of the
  ` vector
  */\
  tmpVU32 =\
     max_vectU32_retVectU32(\
        (inVectU32), \
        srvect_vectU32_retVectU32((inVectU32), 4)\
  );\
  tmpVU32 =\
     max_vectU32_retVectU32(\
        (inVectU32), \
        srvect_vectU32_retVectU32((inVectU32), 2)\
  );\
  (tmpVU32) =\
     max_vectU32_retVectU32(\
        (tmpVU32), \
        srvect_vectU32_retVectU32((tmpVU32), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retU32 = extract_vectU32_U32((tmpVU32), 0);\
   /*extract 32U is a double extract for SSE2 (not SSE4)*/\
}) /*hmaxVectU32*/

/*+++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-06 Cat-08:
+  - Horizontal min (same as hmax, except use min)
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
}) /*hminVectI8*/


#define hmin_vectU8_retU8(retU8, inVectU8) ({\
  vectU8 tmpVU8; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
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
     (uint8_t) \
     extract_vectI16_I16(\
        cast_vectU8_to_VectI16((tmpVU8)),\
        0\
  );\
}) /*hminVectU8*/

#define hmin_vectI16_retI16(retI16, inVectI16) ({\
  vectI16 tmpVI16; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
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
  retI16 = extract_vectI16_I16((tmpVI16), 0);\
}) /*hminVectI16*/

#define hmin_vectU16_retU16(retU16, inVectU16) ({\
  vectU16 tmpVU16; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
  tmpVU16 =\
     min_vectU16_retVectU16(\
        (inVectU16), \
        srvect_vectU16_retVectU16((inVectU16), 4)\
  );\
  (tmpVU16) =\
     min_vectU16_retVectU16(\
        (tmpVU16), \
        srvect_vectU16_retVectU16((tmpVU16), 2)\
  );\
  (tmpVU16) =\
     min_vectU16_retVectU16(\
        (tmpVU16), \
        srvect_vectU16_retVectU16((tmpVU16), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retU16 = extract_vectU16_U16((tmpVU16), 0);\
}) /*minVectU16*/

#define hmin_vectI32_retI32(retI32, inVectI32) ({\
  vectI32 tmpVI32; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
  tmpVI32 =\
     min_vectI32_retVectI32(\
        (inVectI32), \
        srvect_vectI32_retVectI32((inVectI32), 2)\
  );\
  (tmpVI32) =\
     min_vectI32_retVectI32(\
        (tmpVI32), \
        srvect_vectI32_retVectI32((tmpVI32), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retI32 = extract_vectI32_I32((tmpVI32), 0);\
   /*extract 32I is a double extract for SSE2 (not SSE4)*/\
}) /*hminVectI32*/

#define hmin_vectU32_retU32(retU32, inVectU32) ({\
  vectU32 tmpVU32; /*To avoid overwriting input vector*/\
  \
  /*Find and move the min value to the lower part of the
  ` vector
  */\
  tmpVU32 =\
     min_vectU32_retVectU32(\
        (inVectU32), \
        srvect_vectU32_retVectU32((inVectU32), 2)\
  );\
  (tmpVU32) =\
     min_vectU32_retVectU32(\
        (tmpVU32), \
        srvect_vectU32_retVectU32((tmpVU32), 1)\
  );\
  \
  /*Extract the mas value (now at position 0)*/\
  retU32 = extract_vectU32_U32((tmpVU32), 0);\
   /*extract 32U is a double extract for SSE2 (not SSE4)*/\
}) /*hminVectU32*/


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
  ((vectI16) (inVectI8))

#define cast_vectI8_to_VectI32(inVectI8) \
  ((vectI32) (inVectI8))

#define cast_vectI16_to_VectI8(inVectI16) \
  ((vectI8) (inVectI16))

#define cast_vectI16_to_VectI32(inVectI16) \
  ((vectI32) (inVectI16))

#define cast_vectI32_to_VectI8(inVectI32) \
  ((vectI8) (inVectI32))

#define cast_vectI32_to_VectI16(inVectI32) \
  ((vectI16) (inVectI32))

/*Pointer casts*/

#define ptrcast_vectI8_to_VectI16(inVectI8) \
  ((vectI16 *) (inVectI8))

#define ptrcast_vectI8_to_VectI32(inVectI8) \
  ((vectI32 *) (inVectI8))

#define ptrcast_vectI16_to_VectI8(inVectI16) \
  ((vectI8 *) (inVectI16))

#define ptrcast_vectI16_to_VectI32(inVectI16) \
  ((vectI32 *) (inVectI16))

#define ptrcast_vectI32_to_VectI8(inVectI32) \
  ((vectI8 *) (inVectI32))

#define ptrcast_vectI32_to_VectI16(inVectI32) \
  ((vectI16 *) (inVectI32))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-02:
+  - Convert unsigned vectors to unsigned vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectU8_to_VectU16(inVectU8) \
  ((vectU16) (inVectU8))

#define cast_vectU8_to_VectU32(inVectU8) \
  ((vectU32) (inVectU8))

#define cast_vectU16_to_VectU8(inVectU16) \
  ((vectU8) (inVectU16))

#define cast_vectU16_to_VectU32(inVectU16) \
  ((vectU32) (inVectU16))

#define cast_vectU32_to_VectU8(inVectU32) \
  ((vectU8) (inVectU32))

#define cast_vectU32_to_VectU16(inVectU32) \
  ((vectU16) (inVectU32))

/*Pointer casts*/

#define ptrcast_vectU8_to_VectU16(inVectU8) \
  ((vectU16 *) (inVectU8))

#define ptrcast_vectU8_to_VectU32(inVectU8) \
  ((vectU32 *) (inVectU8))

#define ptrcast_vectU16_to_VectU8(inVectU16) \
  ((vectU8 *) (inVectU16))

#define ptrcast_vectU16_to_VectU32(inVectU16) \
  ((vectU32 *) (inVectU16))

#define ptrcast_vectU32_to_VectU8(inVectU32) \
  ((vectU8 *) (inVectU32))

#define ptrcast_vectU32_to_VectU16(inVectU32) \
  ((vectU16 *) (inVectU32))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-03:
+  - Convert unsigned vectors to signed vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectU8_to_VectI8(inVectU8) \
  ((vectI8) (inVectU8))

#define cast_vectU8_to_VectI16(inVectU8) \
  ((vectI16) (inVectU8))

#define cast_vectU8_to_VectI32(inVectU8) \
  ((vectI32) (inVectU8))

#define cast_vectU16_to_VectI8(inVectU16) \
  ((vectI8) (inVectU16))

#define cast_vectU16_to_VectI16(inVectU16) \
  ((vectI16) (inVectU16))

#define cast_vectU16_to_VectI32(inVectU16) \
  ((vectI32) (inVectU16))

#define cast_vectU32_to_VectI8(inVectU32) \
  ((vectI8) (inVectU32))

#define cast_vectU32_to_VectI16(inVectU32) \
  ((vectI16) (inVectU32))

#define cast_vectU32_to_VectI32(inVectU32) \
  ((vectI32) (inVectU32))

/*Pointer casts*/

#define ptrcast_vectU8_to_VectI8(inVectU8) \
  ((vectI8 *) (inVectU8))

#define ptrcast_vectU8_to_VectI16(inVectU8) \
  ((vectI16 *) (inVectU8))

#define ptrcast_vectU8_to_VectI32(inVectU8) \
  ((vectI32 *) (inVectU8))

#define ptrcast_vectU16_to_VectI8(inVectU16) \
  ((vectI8 *) (inVectU16))

#define ptrcast_vectU16_to_VectI16(inVectU16) \
  ((vectI16 *) (inVectU16))

#define ptrcast_vectU16_to_VectI32(inVectU16) \
  ((vectI32 *) (inVectU16))

#define ptrcast_vectU32_to_VectI8(inVectU32) \
  ((vectI8 *) (inVectU32))

#define ptrcast_vectU32_to_VectI16(inVectU32) \
  ((vectI16 *) (inVectU32))

#define ptrcast_vectU32_to_VectI32(inVectU32) \
  ((vectI32 *) (inVectU32))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-04:
+  - Convert signed vectors to unsigned vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define cast_vectI8_to_VectU8(inVectI8) \
  ((vectU8) (inVectI8))

#define cast_vectI8_to_VectU16(inVectI8) \
  ((vectU16) (inVectI8))

#define cast_vectI8_to_VectU32(inVectI8) \
  ((vectU32) (inVectI8))

#define cast_vectI16_to_VectU8(inVectI16) \
  ((vectU8) (inVectI16))

#define cast_vectI16_to_VectU16(inVectI16) \
  ((vectU16) (inVectI16))

#define cast_vectI16_to_VectU32(inVectI16) \
  ((vectU32) (inVectI16))

#define cast_vectI32_to_VectU8(inVectI32) \
  ((vectU8) (inVectI32))

#define cast_vectI32_to_VectU16(inVectI32) \
  ((vectU16) (inVectI32))

#define cast_vectI32_to_VectU32(inVectI32) \
  ((vectU32) (inVectI32))

/*Pointer casts*/

#define ptrcast_vectI8_to_VectU8(inVectI8) \
  ((vectU8 *) (inVectI8))

#define ptrcast_vectI8_to_VectU16(inVectI8) \
  ((vectU16 *) (inVectI8))

#define ptrcast_vectI8_to_VectU32(inVectI8) \
  ((vectU32 *) (inVectI8))

#define ptrcast_vectI16_to_VectU8(inVectI16) \
  ((vectU8 *) (inVectI16))

#define ptrcast_vectI16_to_VectU16(inVectI16) \
  ((vectU16 *) (inVectI16))

#define ptrcast_vectI16_to_VectU32(inVectI16) \
  ((vectU32 *) (inVectI16))

#define ptrcast_vectI32_to_VectU8(inVectI32) \
  ((vectU8 *) (inVectI32))

#define ptrcast_vectI32_to_VectU16(inVectI32) \
  ((vectU16 *) (inVectI32))

#define ptrcast_vectI32_to_VectU32(inVectI32) \
  ((vectU32 *) (inVectI32))

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++\
+ Sec-02 Sub-07 Cat-05:
+  - Convert masks to vectors
\*+++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#define mask8_retVectI8(inMask8) ((vectI8) (inMask8))
#define mask16_retVectI16(inMask16) ((vectI16) (inMask16))
#define mask32_retVectI32(inMask32) ((vectI32) (inMask32))

#define mask8_retVectU8(inMask8) ((vectU8) (inMask8))
#define mask16_retVectU16(inMask16) ((vectU16) (inMask16))
#define mask32_retVectU32(inMask32) ((vectU32) (inMask32))

#endif
