/*
` TODO: NEON Test if vreinterpret can cast pointers
*/

/*#########################################################
# Name: vectorWrap
# Use:
#  - Is a wraper for SSE2 equivlent instruction sets
#    for SSE2, AVX2, and NEON SMID instructions. The idea
#    is to create a series of instructions that the user
#    can depend on being efficent across systems. Commands
#    that might not be efficent are noted in the Functions
#    to watch out for section.
#  - This library uses macros to ensure that the code
#    compiles cleanly every time without the use of -O2 or
#    -O3 flags. That being said, you will still likely use
#    these flags.
#  - SIMDe is likely a better and more robust option.
#  - This code does not compile to scalar, but the
#    vectorBits and vectorBytes are assigned values.
# Compiling for gcc:
#  - sse2: -DSSE2 (Works for all x86 cpus)
#  - sse4: -DSSE4 (better, most x86 systems likely have)
#  - avx2: -DAVX2 -mavx2
#  - avx512: -DAVX512 -mavx2bw (only new intel cpus)
#  - neon: -DNEON -mfpu=neon-vfpv4 or -mfpu=neon (arm)
#  - neon64: -DNEON64 -mfpu=neon-vfpv4 or -mfpu=neon (arm)
#    o I do not think you need this, but it is here just
#      in case.
# Supporting functions:
# Datatype:
#  - vectI8:
#    o Data type for working with chars or int8_t
#  - vectI16:
#    o Data type for working with shorts or int16_t
#  - vectI32:
#    o Data type for working with ints (32 bit is longs) or
#      prefrebly int32_t
#  - vectU8:
#    o Data type for working with unsigned chars or uint8_t
#  - vectU16:
#    o Data type for working with unsigned shorts or
#      uint16_t
#  - vectU32:
#    o Data type for working with unsigned ints (32 bit OS
#      is unsigned longs) or uint32_t
#  - mask8:
#    o Holds result form comparisons (ifs) comparing vectI8
#  - mask16:
#    o Holds result form comparisons (ifs) comparing
#      two vectI16
#  - mask32:
#    o Holds result form comparisons (ifs) comparing
#      two vectI32
# Definitions:
#   - defVectBits:
#     o The number of bits a vector holds.
#     o This goes to 8 when compiled to scalar
#   - defVectBytes:
#     o The number of bytes a vector holds
#     o This goes to 1 when compiled to scalar
#   - defNum8BitElms:
#     o Number of 8 bit elements in a vector.
#     o This goes to 1 when compiling a scalar.
#   - defNum16BitElms:
#     o Number of 16 bit elements in a vector.
#     o This goes to 1 when compiling a scalar.
#   - defNum32BitElms:
#     o Number of 32 bit elements in a vector.
#     o This goes to 1 when compiling a scalar.
#   - defMvMask8
#     o This is the value that is expected to be returned
#        when store_mask8_retUL is used on a mask that
#        has all values true (intel movemask).
#     o For NEON this is a long.
#   - defMvMask16
#     o This is the value that is expected to be returned
#        when store_mask16_retUL is used on a mask that
#        has all values true (intel movemask).
#     o For NEON this is a long.
#   - defMvMask32
#     o This is the value that is expected to be returned
#        when store_mask32_retUL is used on a mask that
#        has all values true (intel movemask).
#     o For NEON this is a long.
# Structures:
#   o vectMemBlock:
#     - Holds buffer for vectors to work with (alnPtr)
#     - This is used to avoid alignment issues for
#       AVX2/AVX512
#     - macros/functions for this structure are in sec-04
#     - How to use:
#       bufferSize = 1024 + defVectBytes;
#       vectMemBlock *vectBuff = initVectMemBlock();
#       mallocVectMemBlock(vectBuff, bufferSize);
#       vectXX vector =
#          load_XXAry_retVectXX(vectBuff->alnPtr);
#       When finshed: freeMemBlock(vectBuff);
# Wrapper functions:
#  Nameing method:
#   o Each function is named by task, input, and output
#   o funciton_input_retOutput
#     - Ary is for array
#     - I or U is for signed or unsigned value
#       o The number after I or U indicates the number of
#         bit should be in the input number.
#     - vect, Vect, or V is for vector
#     - mask or Mask is for a mask
#  Functions to watch out for:
#   - There are a few functions that are not one to one
#     functions. Some of these are easy to get around and
#     still worth doing, while others, might not be worth
#     using at all.
#    - List of problematic functions:
#      o unsigned less than or greater than comparisons
#        for SSE2 (5 extra opterations + 3 additional
#        vectors).
#        - This is limited to SSE2 only.
#      o unsigned 16 and 32 bit maximums and minimums for
#        SSE2. 8 operations + 3 (likely not 4?)
#        extra vectors.
#        - This is limited to SSE2 only.
#        - This is worth it for unsigned 16bit, but is
#          a bit borderline for 32bit
#          o You will need to load in the input vectors and
#            store the maximum result, this makes this
#            function a bit borderline in costs.
#    - List of functions requiring extra steps
#      o NEON does not have a movemask equivlent and AVX512
#        does not return a vector for comparisions. To get
#        around this I added in something close to a
#        movemask command for NEON. Then I added a
#        fix_maskX_popcount, which will convert the
#        population counts for SSE, AVX2, AVX512, and NEON.
#      o workflow:
#        - cmpxx_vectXX_retMaskX(x, y);
#          xUL = store_maskX_retUL;
#          xUL = fix_maskX_popcount(xUL);
#        - NEON move mask command from:
#          https://community.arm.com/arm-community-blogs/b/infrastructure-solutions-blog/posts/porting-x86-vector-bitmask-optimizations-to-arm-neon
#      o Basicly, be carefull when manipulating masks. They
#        are not consistent.
#  Functions:
#  - Functions either take ro return a minVectIX or vactIX.
#  - Any parameter that is named ret is modified. Other
#    than those functions all functions return a value.
#  - The sec-xx sub-xx cat-xx refers to the functions
#    location in vectWrapAVX512.h, vectWrapAVX2.h,
#    vectWrapSSE.h, vectWrapNEON.h, or vectWrapNEON64.h.
#  - Input (sec-02 sub-01)
#    o Load in aligned data (sec-02 sub-01 cat-01)
#      - load_I8Ary_retVectI8(arrayI8):
#        o Loads in an array of int8_t into a vector
#      - load_I16Ary_retVectI16(arrayI16):
#        o Loads in an array of int16_t into a vector
#      - load_I32Ary_retVectI32(arrayI32):
#        o Loads in an array of int32_t into a vector
#      - load_U8Ary_retVectU8(arrayU8):
#        o Loads in an array of uint8_t into a vector
#      - load_U16Ary_retVectU16(arrayU16):
#        o Loads in an array of uint16_t into a vector
#      - load_U32Ary_retVectU32(arrayU32):
#        o Loads in an array of uint32_t into a vector
#    o Load in unaligned data (sec-02 sub-01 cat-01)
#      - loadu_I8Ary_retVectI8(arrayI8):
#        o Loads in an array of int8_t into a vector
#      - loadu_I16Ary_retVectI16(arrayI16):
#        o Loads in an array of int16_t into a vector
#      - loadu_I32Ary_retVectI32(arrayI32):
#        o Loads in an array of int32_t into a vector
#      - loadu_U8Ary_retVectU8(arrayU8):
#        o Loads in an array of uint8_t into a vector
#      - loadu_U16Ary_retVectU16(arrayU16):
#        o Loads in an array of uint16_t into a vector
#      - loadu_U32Ary_retVectU32(arrayU32):
#        o Loads in an array of uint32_t into a vector
#    o Make a vector of zeros (sec-02 sub-01 cat-03)
#      - NOTE: These may be slower on NEON. I have to use
#        vdup_n_s8(0) and I am not sure the cost of this.
#      - zero_retVectI8():
#        o Returns a vectI8 full of zeros
#      - zero_retVectI16():
#        o Returns a vectI16 full of zeros
#      - zero_retVectI32():
#        o Returns a vectI32 full of zeros
#      - zero_retVectU8():
#        o Returns a vectU8 full of zeros
#      - zero_retVectU16():
#        o Returns a vectU16 full of zeros
#      - zero_retVectU32():
#        o Returns a vectU32 full of zeros
#    o Make vector from single value (sec-02 sub-01 cat-04)
#      - set1_I8_retVectI8(valI8);
#        o Returns an vectI8 with every element being valI8
#      - set1_I16_retVectI16(valI16);
#        o Returns an vectI16 with every element being
#          valI16
#      - set1_I32_retVectI32(valI32);
#        o Returns an vectI32 with every element being
#          valI32
#      - set1_U8_retVectU8(valU8);
#        o Returns an vectU8 with every element being valU8
#      - set1_U16_retVectU16(valU16);
#        o Returns an vectU16 with every element being
#          valU16
#      - set1_U32_retVectU32(valU32);
#        o Returns an vectU32 with every element being
#          valU32
#    o Insert element into an vector (sec-02 sub-01 cat-04)
#      - insert_I16_retVectI16(vectI16,insertI16,position);
#        o Inserts insertI16 (int16_t) into vectI16 at
#          position.
#      - insert_I32_retVectI32(vectI32,insertI32,position);
#        o Inserts insertI32 (int32_t) into vectI32 at
#          position.
#        o WARNING: ON SSE2, but not SSE4 this does
#          insert_I16_retVectI16 twice. So, this is slower
#          on SSE2 (latency of 4 instead of 2).
#      - insert_U16_retVectU16(vectU16,insertU16,position);
#        o Unserts insertU16 (uint16_t) into vectU16 at
#          position.
#      - insert_U32_retVectU32(vectU32,insertU32,position);
#        o Unserts insertU32 (uint32_t) into vectU32 at
#          position.
#        o WARNING: ON SSE2, but not SSE4 this does
#          insert_I16_retVectI16 twice and needs an 
#          additional bit shift. So, this is slower on
#          SSE2 (latency of 4 instead of 2).
#  - Output (sec-02 sub-02)
#    o Put vector into an array (sec-02 sub-02 cat-01)
#      - store_vectI8_retAryI8(arrayI8, inVectI8)
#        o Stores contents of inVectI8 into an array of
#          int8_t's
#      - store_vectI16_retAryI16(arrayI16, inVectI16)
#        o Stores contents of inVectI16 into an array of
#          int16_t's
#      - store_vectI32_retAryI32(arrayI32, inVectI32)
#        o Stores contents of inVectI32 into an array of
#          int32_t's
#      - store_vectU8_retAryU8(arrayU8, inVectU8)
#        o Stores contents of inVectU8 uinto an array of
#          uint8_t's
#      - store_vectU16_retAryU16(arrayU16, inVectU16)
#        o Stores contents of inVectU16 uinto an array of
#          uint16_t's
#      - store_vectU32_retAryU32(arrayU32, inVectU32)
#        o Stores contents of inVectU32 uinto an array of
#          uint32_t's
#    o Put vector in unaligned array (sec-02 sub-02 cat-02)
#      - storeu_vectI8_retAryI8(arrayI8, inVectI8)
#        o Stores contents of inVectI8 into an array of
#          int8_t's (does not need to be aligned)
#      - storeu_vectI16_retAryI16(arrayI16, inVectI16)
#        o Stores contents of inVectI16 into an array of
#          int16_t's (does not need to be aligned)
#      - storeu_vectI32_retAryI32(arrayI32, inVectI32)
#        o Stores contents of inVectI32 into an array of
#          int32_t's (does not need to be aligned)
#      - storeu_vectU8_retAryU8(arrayU8, inVectU8)
#        o Stores contents of inVectU8 uinto an array of
#          uint8_t's (does not need to be aligned)
#      - storeu_vectU16_retAryU16(arrayU16, inVectU16)
#        o Stores contents of inVectU16 uinto an array of
#          uint16_t's (does not need to be aligned)
#      - storeu_vectU32_retAryU32(arrayU32, inVectU32)
#        o Stores contents of inVectU32 uinto an array of
#          uint32_t's (does not need to be aligned)
#    o Store masks into longs (sec-02 sub-02 cat-03)
#      - store_mask8_retUL(returnULong, inMask8)
#        o Converts an mask8 into an unsinged long
#      - store_mask16_retUL(returnULong, inMask16)
#        o Converts an mask16 into an unsinged long
#      - store_mask32_retUL(returnULong, inMask32)
#        o Converts an mask32 into an unsinged long
#    o Extract an element (sec-02 sub-02 cat-04)
#      - extract_vectI8_I8(inVectI8, posI)\
#        o Extracts the 8 bit element at posI from a vector
#        o Note for AVX512 this reads in posI, shuffles,
#          and then extracts (latency is over 8)
#        o Note: for SSE2 this does a 16 bit extract and
#          then shifts/casts out the incorrect bits (4 op)
#      - extract_vectU8_U8(inVectU8, posI)\
#        o Extracts the 8 bit element at posI from a vector
#        o Note for AVX512 this reads in posI, shuffles,
#          and then extracts (latency is over 8)
#        o Note: for SSE2 this does a 16 bit extract and
#          then shifts/casts out the incorrect bits (4 op)
#      - extract_vectI16_I16(inVectI16, posI)\
#        o Extracts 16 bit element at posI from a vector
#        o Note for AVX512 this reads in posI, shuffles,
#          and then extracts (latency is over 8)
#      - extract_vectU16_U16(inVectU16, posI)\
#        o Extracts 16 bit element at posI from a vector
#        o Note for AVX512 this reads in posI, shuffles,
#          and then extracts (latency is over 8)
#      - extract_vectI32_I32(inVectI32, posI)\
#        o Note for AVX512 this reads in posI, shuffles,
#          and then extracts (latency is over 8)
#        o Note for SSE2 this does a double extract
#      - extract_vectU32_U32(inVectU32, posI)\
#        o Note for AVX512 this reads in posI, shuffles,
#          and then extracts (latency is over 8)
#        o Note for SSE2 this does a double extract
#  - Comparisons (sec-02 sub-03)
#    o signed comparisons (sec-02 sub-03 cat-01)
#      - cmpeq_vectI8_retMask8(xVectI8, yVectI8)
#        o Returns an mask8 holding -1 were xVectI8 is
#          equal to yVectI8 or 0's were not equal.
#      - cmplt_vectI8_retMask8(xVectI8, yVectI8)
#        o Returns an mask8 holding -1 were xVectI8 is
#          less than yVectI8 or 0's were equal or >.
#      - cmpgt_vectI8_retMask8(xVectI8, yVectI8)
#        o Returns an mask8 holding -1 were xVectI8 is
#          greater than yVectI8 or 0's were equal or <.
#      - cmpeq_vectI16_retMask16(xVectI16, yVectI16)
#        o Returns an mask16 holding -1 were xVectI16 is
#          equal to yVectI16 or 0's were not equal.
#      - cmplt_vectI16_retMask16(xVectI16, yVectI16)
#        o Returns an mask16 holding -1 were xVectI16 is
#          less than yVectI16 or 0's were equal or >.
#      - cmpgt_vectI16_retMask16(xVectI16, yVectI16)
#        o Returns an mask16 holding -1 were xVectI16 is
#          greater than yVectI16 or 0's were equal or <.
#      - cmpeq_vectI32_retMask32(xVectI32, yVectI32)
#        o Returns an mask32 holding -1 were xVectI32 is
#          equal to yVectI32 or 0's were not equal.
#      - cmplt_vectI32_retMask32(xVectI32, yVectI32)
#        o Returns an mask32 holding -1 were xVectI32 is
#          less than yVectI32 or 0's were equal or >.
#      - cmpgt_vectI32_retMask32(xVectI32, yVectI32)
#        o Returns an mask32 holding -1 were xVectI32 is
#          greater than yVectI32 or 0's were equal or <.
#    o unsinged equal to comparisons (sec-02 sub-03 cat-02)
#      - cmpeq_vectU8_retMask8(xVectU8, yVectU8)
#        o Returns an mask8 holding -1 were xVectU8 is
#          equal to yVectU8 or 0's were not equal.
#      - cmpeq_vectU16_retMask16(xVectU16, yVectU16)
#        o Returns an mask16 holding -1 were xVectU16 is
#          equal to yVectU16 or 0's were not equal.
#      - cmpeq_vectU32_retMask32(xVectU32, yVectU32)
#        o Returns an mask32 holding -1 were xVectU32 is
#          equal to yVectU32 or 0's were not equal.
#    o unsinged less than comparison (sec-02 sub-03 cat-03)
#      - WARNING: Thes functions take 9 opertations and two
#        additional vectors in SSE2 (not SSE4). These
#        functions should be avoided when possible.
#      - cmplt_vectU8_retMask8(xVectU8, yVectU8)
#        o Returns an mask8 holding -1 were xVectU8 is
#          less than yVectU8 or 0's were not less than.
#      - cmplt_vectU16_retMask16(xVectU16, yVectU16)
#        o Returns an mask16 holding -1 were xVectU16 is
#          less than yVectU16 or 0's were not less than.
#      - cmplt_vectU32_retMask32(xVectU32, yVectU32)
#        o Returns an mask32 holding -1 were xVectU32 is
#          less than to yVectU32 or 0's were not less than.
#    o unsinged > comparison (sec-02 sub-03 cat-04)
#      - WARNING: These functions take 9 opertations & two
#        additional vectors in SSE2 (not SSE4). These
#        functions should be avoided when possible.
#      - cmpgt_vectU8_retMask8(xVectU8, yVectU8)
#        o Returns an mask8 holding -1 were xVectU8 is
#          > yVectU8 or 0's were not >.
#      - cmpgt_vectU16_retMask16(xVectU16, yVectU16)
#        o Returns an mask16 holding -1 were xVectU16 is
#          > yVectU16 or 0's were not >.
#      - cmpgt_vectU32_retMask32(xVectU32, yVectU32)
#        o Returns an mask32 holding -1 were xVectU32 is
#          > to yVectU32 or 0's were not >.
#    o Fix differences in pop counts (sec-02 sub-03 cat-05)
#      - fix_mask8_popcount(inUL)
#        o Adjust for doing a population count with an
#          mask8
#      - fix_mask16_popcount(inUL)
#        o Adjust for doing a population count with an
#          mask16
#      - fix_mask32_popcount(inUL)
#        o Adjust for doing a population count with an
#          mask32
#  - Logical functions (sec-02 sub-04)
#    o andNot functions (sec-02 sub-04 cat-01)
#      - andNot_vectI8_retVectI8(notVectI8, inVectI8)
#        o do an and operations on the not(notVectI8) and
#          inVectI8
#      - andNot_vectI16_retVectI16(notVectI16, inVectI16)
#        o do an and operations on the not(notVectI16) and
#          inVectI16
#      - andNot_vectI32_retVectI32(notVectI32, inVectI32)
#        o do an and operations on the not(notVectI32) and
#          inVectI32
#      - andNot_vectU8_retVectU8(notVectU8, inVectU8)
#        o do an and operations on the not(notVectU8) and
#          inVectU8
#      - andNot_vectU16_retVectU16(notVectU16, inVectU16)
#        o do an and operations on the not(notVectU16) and
#          inVectU16
#      - andNot_vectU32_retVectU32(notVectU32, inVectU32)
#        o do an and operations on the not(notVectU32) and
#          inVectU32
#    o and functions (sec-02 sub-04 cat-02)
#      - and_vectI8_retVectI8(xVectI8, yVectI8)
#        o do an and operation on the xVectI8 and yVectI8
#      - and_vectI16_retVectI16(xVectI16, yVectI16)
#        o do an and operation on xVectI16 and yVectI16
#      - and_vectI32_retVectI32(xVectI32, yVectI32)
#        o do an and operation on xVectI32 and yVectI32
#      - and_vectU8_retVectU8(xVectU8, yVectU8)
#        o do an and operation on xVectU8 and yVectU8
#      - and_vectU16_retVectU16(xVectU16, yVectU16)
#        o do an and operation on xVectU16 and yVectU16
#      - and_vectU32_retVectU32(xVectU32, yVectU32)
#        o do an and operation on xVectU32 and yVectU32
#    o or functions (sec-02 sub-04 cat-03)
#      - or_vectI8_retVectI8(xVectI8, yVectI8)
#        o do an or operation on the xVectI8 or yVectI8
#      - or_vectI16_retVectI16(xVectI16, yVectI16)
#        o do an or operation on xVectI16 or yVectI16
#      - or_vectI32_retVectI32(xVectI32, yVectI32)
#        o do an or operation on xVectI32 or yVectI32
#      - or_vectU8_retVectU8(xVectU8, yVectU8)
#        o do an or operation on xVectU8 or yVectU8
#      - or_vectU16_retVectU16(xVectU16, yVectU16)
#        o do an or operation on xVectU16 or yVectU16
#      - or_vectU32_retVectU32(xVectU32, yVectU32)
#        o do an or operation on xVectU32 and yVectU32
#    o xor functions (sec-02 sub-04 cat-04)
#      - xor_vectI8_retVectI8(xVectI8, yVectI8)
#        o do an xor operation on the xVectI8 xor yVectI8
#      - xor_vectI16_retVectI16(xVectI16, yVectI16)
#        o do an xor operation on xVectI16 xor yVectI16
#      - xor_vectI32_retVectI32(xVectI32, yVectI32)
#        o do an xor operation on xVectI32 xor yVectI32
#      - xor_vectU8_retVectU8(xVectU8, yVectU8)
#        o do an xor operation on xVectU8 xor yVectU8
#      - xor_vectU16_retVectU16(xVectU16, yVectU16)
#        o do an xor operation on xVectU16 xor yVectU16
#      - xor_vectU32_retVectU32(xVectU32, yVectU32)
#        o do an xor operation on xVectU32 and yVectU32
#  - Bit shifts (sec-02 sub-05)
#    o right shifts (sec-02 sub-05 cat-01)
#      - sr_vectI16_retVectI16(inVectI16, numBitsI)
#        o Shift each int16_t element in inVectI16 right by
#          numBitsI.
#      - sr_vectI32_retVectI32(inVectI32, numBitsI)
#        o Shift each int32_t element in inVectI32 right by
#          numBitsI.
#      - sr_vectU16_retVectU16(inVectU16, numBitsU)
#        o Shift each uint16_t element in inVectU16 right
#          by numBitsU.
#      - sr_vectU32_retVectU32(inVectU32, numBitsU)
#        o Shift each uint32_t element in inVectU32 right
#          by numBitsU.
#    o right shift (entire vector) (sec-02 sub-05 cat-02)
#      - srvect_vectI8_retVectI8(inVectI8, numBytesI)
#        o Shift inVectI8 right by numBytesI bytes
#      - srvect_vectI16_retVectI16(inVectI16, numShortsI)
#        o Shift inVectI16 right by numShortsI shorts
#      - srvect_vectI32_retVectI32(inVectI32, numIntsI)
#        o Shift inVectI32 right by numIntsI ints
#      - srvect_vectU8_retVectU8(inVectU8, numBytesU)
#        o Shift inVectU8 right by numBytesU bytes
#      - srvect_vectU16_retVectU16(inVectU16, numShortsU)
#        o Shift inVectU16 right by numShortsU shorts
#      - srvect_vectU32_retVectU32(inVectU32, numUntsU)
#        o Shift inVectU32 right by numUntsU ints
#    o left shifts (sec-02 sub-05 cat-03)
#      - sr_vectI16_retVectI16(inVectI16, numBitsI)
#        o Shift each int16_t element in inVectI16 left by
#          numBitsI.
#      - sr_vectI32_retVectI32(inVectI32, numBitsI)
#        o Shift each int32_t element in inVectI32 left by
#          numBitsI.
#      - sr_vectU16_retVectU16(inVectU16, numBitsU)
#        o Shift each uint16_t element in inVectU16 left
#          by numBitsU.
#      - sr_vectU32_retVectU32(inVectU32, numBitsU)
#        o Shift each uint32_t element in inVectU32 left
#          by numBitsU.
#    o left shift (entire vector) (sec-02 sub-05 cat-04)
#      - slvect_vectI8_retVectI8(inVectI8, numBytesI)
#        o Shift inVectI8 left by numBytesI bytes
#      - slvect_vectI16_retVectI16(inVectI16, numShortsI)
#        o Shift inVectI16 left by numShortsI shorts
#      - slvect_vectI32_retVectI32(inVectI32, numIntsI)
#        o Shift inVectI32 left by numIntsI ints
#      - slvect_vectU8_retVectU8(inVectU8, numBytesU)
#        o Shift inVectU8 left by numBytesU bytes
#      - slvect_vectU16_retVectU16(inVectU16, numShortsU)
#        o Shift inVectU16 left by numShortsU shorts
#      - slvect_vectU32_retVectU32(inVectU32, numUntsU)
#        o Shift inVectU32 left by numUntsU ints
#  - Math (sec-02 sub-06)
#    o additions (sec-02 sub-06 cat-01)
#      - add_vectI8_retVectI8(xVectI8, yVectI8)
#        o Add together xVectI8 and yVectI8.
#      - add_vectI16_retVectI16(xVectI16, yVectI16)
#        o Add together xVectI16 and yVectI16.
#      - add_vectI32_retVectI32(xVectI32, yVectI32)
#        o Add together xVectI32 and yVectI32.
#    o Saturation additions (sec-02 sub-06 cat-02)
#      - addSat_vectU8_retVectU8(xVectU8, yVectU8)
#        o Do a stauration additions with xVectU8 and
#          yVectU8
#      - addSat_vectU16_retVectU16(xVectU16, yVectU16)
#        o Do a stauration addition with xVectU16 and
#          yVectI16
#      - addSat_vectI8_retVectI8(xVectI8, yVectI8)
#        o Do a stauration additions with xVectI8 and
#          yVectI8
#      - addSat_vectI16_retVectI16(xVectI16, yVectI16)
#        o Do a stauration addition with xVectI16 and
#          yVectI16
#    o subtrations (sec-02 sub-06 cat-03)
#      - sub_vectI8_retVectI8(xVectI8, yVectI8)
#        o Subract xVectI8 from yVectI8.
#      - sub_vectI16_retVectI16(xVectI16, yVectI16)
#        o subrract xVectI16 from yVectI16.
#      - sub_vectI32_retVectI32(xVectI32, yVectI32)
#        o subtract xVectI32 from yVectI32.
#    o Saturation subtrations (sec-02 sub-06 cat-04)
#      - subSat_vectU8_retVectU8(xVectU8, yVectU8)
#        o Do a stauration subtraction with xVectU8 and
#          yVectU8
#      - subSat_vectU16_retVectU16(xVectU16, yVectU16)
#        o Do a stauration subtraction with xVectU16 and
#          yVectU16
#      - subSat_vectI8_retVectI8(xVectI8, yVectI8)
#        o Do a stauration subtraction with xVectI8 and
#          yVectI8
#      - subSat_vectI16_retVectI16(xVectI16, yVectI16)
#        o Do a stauration subtraction with xVectI16 and
#          yVectI16
#    o Max (sec-02 sub-06 cat-05)
#      - max_vectI8_retVectI8(xVectI8, yVectI8)
#        o Find the maximum values between xVect and yVect
#        o Note: For SSE2 (not SSE4) this requires two
#          extra vectors and four operations.
#      - max_vectI16_retVectI16(xVectI16, yVectI16)
#        o Find the maximum values between xVect and yVect
#      - max_vectI32_retVectI32(xVectI32, yVectI32)
#        o Find the maximum values between xVect and yVect
#        o Note: For SSE2 (not SSE4) this requires two
#          extra vectors and four operations.
#      - max_vectU8_retVectU8(xVectU8, yVectU8)
#        o Find the maximum values between xVect and yVect
#      - max_vectU16_retVectU16(xVectU16, yVectU16)
#        o Find the maximum values between xVect and yVect
#        o WARINING: For SSE2 (not SSE4) this requires 
#          11 operations and three extra vectors
#      - max_vectU32_retVectU32(xVectU32, yVectU32)
#        o Find the maximum values between xVect and yVect
#        o WARINING: For SSE2 (not SSE4) this requires 
#          11 operations and three extra vectors.
#    o Min (sec-02 sub-06 cat-06)
#      - min_vectI8_retVectI8(xVectI8, yVectI8)
#        o Find the minimum values between xVect and yVect
#        o Note: For SSE2 (not SSE4) this requires two
#          extra vectors and four operations.
#      - min_vectI16_retVectI16(xVectI16, yVectI16)
#        o Find the minimum values between xVect and yVect
#      - min_vectI32_retVectI32(xVectI32, yVectI32)
#        o Find the minimum values between xVect and yVect
#        o Note: For SSE2 (not SSE4) this requires two
#          extra vectors and four operations.
#      - min_vectU8_retVectU8(xVectU8, yVectU8)
#        o Find the minimum values between xVect and yVect
#      - min_vectU16_retVectU16(xVectU16, yVectU16)
#        o Find the minimum values between xVect and yVect
#        o WARINING: For SSE2 (not SSE4) this requires 
#          12 operations and two extra vectors
#      - min_vectU32_retVectU32(xVectU32, yVectU32)
#        o Find the minimum values between xVect and yVect
#        o WARINING: For SSE2 (not SSE4) this requires 
#          12 operations and two extra vectors.
#    o Horizontal max (max in vect) (sec-02 sub-06 cat-07)
#      - hmax_vectI8_retI8(retI8, vectI8)
#        o Modifies retI8 to hold max value in vectI8
#        o Calls 8 opterions (+1 AVX2, +2 AVX512) +
#          1 extract. It also uses a temporary vector.
#      - hmax_vectU8_retU8(retU8, vectU8)
#        o Modifies retI8 to hold max value in vectI8
#        o Calls 9 opterions (+1 AVX2, +2 AVX512) +
#          1 extract. It also uses a temporary vector.
#      - hmax_vectI16_retI16(retI16, vectI16)
#        o Modifies retI16 to hold max value in vectI16
#        o Calls 7 opterions (+1 AVX2, +2 AVX512) +
#          1 extract. It also uses a temporary vector.
#      - hmax_vectU16_retU16(retU16, vectU16)
#        o Modifies retI16 to hold max value in vectI16
#        o Calls 7 opterions (+1 AVX2, +2 AVX512) +
#          1 extract. It also uses a temporary vector.
#      - hmax_vectI32_retI32(retI32, vectI32)
#        o Modifies retI32 to hold max value in vectI32
#        o Calls 5 opterions (+1 AVX2, +2 AVX512?) +
#          1 extract (two extracts for SSE2). It also
#          uses a temporary vector.
#      - hmax_vectU32_retU32(retU32, vectU32)
#        o Modifies retI32 to hold max value in vectI32
#        o Calls 5 opterions (+1 AVX2, +2 AVX512?) +
#          1 extract (two extracts for SSE2). It also
#          uses a temporary vector.
#    o Horizontal min (min in vect) (sec-02 sub-06 cat-07)
#      - hmin_vectI8_retI8(retI8, vectI8)
#        o Modifies retI8 to hold min value in vectI8
#        o Calls 8 opterions (+1 AVX2, +2 AVX512) +
#          1 extract. It also uses a temporary vector.
#      - hmin_vectU8_retU8(retU8, vectU8)
#        o Modifies retI8 to hold min value in vectI8
#        o Calls 9 opterions (+1 AVX2, +2 AVX512) +
#          1 extract. It also uses a temporary vector.
#      - hmin_vectI16_retI16(retI16, vectI16)
#        o Modifies retI16 to hold min value in vectI16
#        o Calls 7 opterions (+1 AVX2, +2 AVX512) +
#          1 extract. It also uses a temporary vector.
#      - hmin_vectU16_retU16(retU16, vectU16)
#        o Modifies retI16 to hold min value in vectI16
#        o Calls 7 opterions (+1 AVX2, +2 AVX512) +
#          1 extract. It also uses a temporary vector.
#      - hmin_vectI32_retI32(retI32, vectI32)
#        o Modifies retI32 to hold min value in vectI32
#        o Calls 5 opterions (+1 AVX2, +2 AVX512?) +
#          1 extract (two extracts for SSE2). It also
#          uses a temporary vector.
#      - hmin_vectU32_retU32(retU32, vectU32)
#        o Modifies retI32 to hold min value in vectI32
#        o Calls 5 opterions (+1 AVX2, +2 AVX512?) +
#          1 extract (two extracts for SSE2). It also
#          uses a temporary vector.
#  - Conversion (sec-02 sub-07)
#    o Signed vectors to signed vect (sec-02 sub-07 cat-01)
#      - These do not compile, but are here to avoid
#        compiler errors (NEON uses these types).
#      - cast_vectI8_to_VectI16(inVectI8)
#        o Convert a vectI8 to a vectI16
#      - cast_vectI8_to_VectI32(inVectI8)
#        o Convert a vectI8 to a vectI32
#      - cast_vectI16_to_VectI8(inVectI16)
#        o Convert a vectI16 to a vectI8
#      - cast_vectI16_to_VectI32(inVectI16)
#        o Convert a vectI16 to a vectI32
#      - cast_vectI32_to_VectI8(inVectI32)
#        o Convert a vectI32 to a vectI8
#      - cast_vectI32_to_VectI16(inVectI32)
#        o Convert a vectI32 to a vectI16
#      - ptrcast_vectI8_to_VectI16(inVectI8)
#        o Convert a vectI8 pionter to a vectI16 pointer
#      - ptrcast_vectI8_to_VectI32(inVectI8)
#        o Convert a vectI8 pointer to a vectI32 pointer
#      - ptrcast_vectI16_to_VectI8(inVectI16)
#        o Convert a vectI16 pointer to a vectI8 pointer
#      - ptrcast_vectI16_to_VectI32(inVectI16)
#        o Convert a vectI16 pointer to a vectI32 pointer
#      - ptrcast_vectI32_to_VectI8(inVectI32)
#        o Convert a vectI32 pointer to a vectI8 pointer
#      - ptrcast_vectI32_to_VectI16(inVectI32)
#        o Convert a vectI32 pionter to a vectI16 pointer
#    o unsigned vect to unsigned vect(sec-02 sub-07 cat-02)
#      - These do not compile, but are here to avoid
#        compiler errors (NEON uses these types).
#      - cast_vectU8_to_VectU16(inVectU8)
#        o Convert a vectU8 to a vectU16
#      - cast_vectU8_to_VectU32(inVectU8)
#        o Convert a vectU8 to a vectU32
#      - cast_vectU16_to_VectU8(inVectU16)
#        o Convert a vectU16 to a vectU8
#      - cast_vectU16_to_VectU32(inVectU16)
#        o Convert a vectU16 to a vectU32
#      - cast_vectU32_to_VectU8(inVectU32)
#        o Convert a vectU32 to a vectU8
#      - cast_vectU32_to_VectU16(inVectU32)
#        o Convert a vectU32 to a vectU16
#      - ptrcast_vectU8_to_VectU16(inVectU8)
#        o Convert a vectU8 pionter to a vectU16 pointer
#      - ptrcast_vectU8_to_VectU32(inVectU8)
#        o Convert a vectU8 pointer to a vectU32 pointer
#      - ptrcast_vectU16_to_VectU8(inVectU16)
#        o Convert a vectU16 pointer to a vectU8 pointer
#      - ptrcast_vectU16_to_VectU32(inVectU16)
#        o Convert a vectU16 pointer to a vectU32 pointer
#      - ptrcast_vectU32_to_VectU8(inVectU32)
#        o Convert a vectU32 pointer to a vectU8 pointer
#      - ptrcast_vectU32_to_VectU16(inVectU32)
#        o Convert a vectU32 pionter to a vectU16 pointer
#    o unsigned vect to signed vect (sec-02 sub-07 cat-03)
#      - These do not compile, but are here to avoid
#        compiler errors (NEON uses these types).
#      - cast_vectU8_to_VectI8(inVectI8)
#        o Convert a vectU8 to a vectI8
#      - cast_vectU8_to_VectI16(inVectI8)
#        o Convert a vectU8 to a vectI16
#      - cast_vectU8_to_VectI32(inVectI8)
#        o Convert a vectU8 to a vectI32
#      - cast_vectU16_to_VectI8(inVectI16)
#        o Convert a vectU16 to a vectI8
#      - cast_vectU16_to_VectI16(inVectI16)
#        o Convert a vectU16 to a vectI16
#      - cast_vectU16_to_VectI32(inVectI16)
#        o Convert a vectU16 to a vectI32
#      - cast_vectU32_to_VectI8(inVectI32)
#        o Convert a vectU32 to a vectI8
#      - cast_vectU32_to_VectI16(inVectI32)
#        o Convert a vectU32 to a vectI16
#      - cast_vectU32_to_VectI32(inVectI32)
#        o Convert a vectU32 to a vectI32
#      - ptrcast_vectU8 pointer_to_VectI8(inVectI8)
#        o Convert a vectU8 pointer to a vectI8 pointer
#      - ptrcast_vectU8 pointer_to_VectI16(inVectI8)
#        o Convert a vectU8 pointer to a vectI16 pointer
#      - ptrcast_vectU8 pointer_to_VectI32(inVectI8)
#        o Convert a vectU8 pointer to a vectI32 pointer
#      - ptrcast_vectU16 pointer_to_VectI8(inVectI16)
#        o Convert a vectU16 pointer to a vectI8 pointer
#      - ptrcast_vectU16 pointer_to_VectI16(inVectI16)
#        o Convert a vectU16 pointer to a vectI16 pointer
#      - ptrcast_vectU16 pointer_to_VectI32(inVectI16)
#        o Convert a vectU16 pointer to a vectI32 pointer
#      - ptrcast_vectU32 pointer_to_VectI8(inVectI32)
#        o Convert a vectU32 pointer to a vectI8 pointer
#      - ptrcast_vectU32 pointer_to_VectI16(inVectI32)
#        o Convert a vectU32 pointer to a vectI16 pointer
#      - ptrcast_vectU32 pointer_to_VectI32(inVectI32)
#        o Convert a vectU32 pointer to a vectI32 pointer
#    o signed vect to unsigned vect (sec-02 sub-07 cat-04)
#      - These do not compile, but are here to avoid
#        compiler errors (NEON uses these types).
#      - cast_vectI8_to_VectU8(inVectU8)
#        o Convert a vectI8 to a vectU8
#      - cast_vectI8_to_VectU16(inVectU8)
#        o Convert a vectI8 to a vectU16
#      - cast_vectI8_to_VectU32(inVectU8)
#        o Convert a vectI8 to a vectU32
#      - cast_vectI16_to_VectU8(inVectU16)
#        o Convert a vectI16 to a vectU8
#      - cast_vectI16_to_VectU16(inVectU16)
#        o Convert a vectI16 to a vectU16
#      - cast_vectI16_to_VectU32(inVectU16)
#        o Convert a vectI16 to a vectU32
#      - cast_vectI32_to_VectU8(inVectU32)
#        o Convert a vectI32 to a vectU8
#      - cast_vectI32_to_VectU16(inVectU32)
#        o Convert a vectI32 to a vectU16
#      - cast_vectI32_to_VectU32(inVectU32)
#        o Convert a vectI32 to a vectU32
#      - ptrcast_vectI8 pointer_to_VectU8(inVectU8)
#        o Convert a vectI8 pointer to a vectU8 pointer
#      - ptrcast_vectI8 pointer_to_VectU16(inVectU8)
#        o Convert a vectI8 pointer to a vectU16 pointer
#      - ptrcast_vectI8 pointer_to_VectU32(inVectU8)
#        o Convert a vectI8 pointer to a vectU32 pointer
#      - ptrcast_vectI16 pointer_to_VectU8(inVectU16)
#        o Convert a vectI16 pointer to a vectU8 pointer
#      - ptrcast_vectI16 pointer_to_VectU16(inVectU16)
#        o Convert a vectI16 pointer to a vectU16 pointer
#      - ptrcast_vectI16 pointer_to_VectU32(inVectU16)
#        o Convert a vectI16 pointer to a vectU32 pointer
#      - ptrcast_vectI32 pointer_to_VectU8(inVectU32)
#        o Convert a vectI32 pointer to a vectU8 pointer
#      - ptrcast_vectI32 pointer_to_VectU16(inVectU32)
#        o Convert a vectI32 pointer to a vectU16 pointer
#      - ptrcast_vectI32 pointer_to_VectU32(inVectU32)
#        o Convert a vectI32 pointer to a vectU32 pointer
#    o Convert masks to vectors (sec-02 sub-07 cat-05)
#      - NOTE: This does have an latency of 2 on AVX512,
#        but is not compiled for AVX2, SSE, or NEON.
#      - mask8_retVectI8(inMask8)
#        o Converts mask8 to a vectI8
#      - mask16_retVectI16(inMask16)
#        o Converts mask16 to a vectI16
#      - mask32_retVectI32(inMask32)
#        o Converts mask32 to a vectI32
#      - mask8_retVectU8(inMask8)
#        o Converts mask8 to a vectU8
#      - mask16_retVectU16(inMask16)
#        o Converts mask16 to a vectU16
#      - mask32_retVectU32(inMask32)
#        o Converts mask32 to a vectU32
# Libraries:
#  o Compiler flag determined
#    i) or "vectWrapAVX512.h" 
#    ii) or "vectWrapAVX2.h"
#    iii) or "vectWrapSSE.h"
#    iv) or "vectWrapNEON.h"
#    v) or "vectWrapNEON64.h"
# C Standard Libraries:
#   - <limits.h>
#   - <stdlib.h>
#   o <stdint.h>
#   o Compiler flaged determined
#     i) or <immintrin.h>
#     ii) or <arm_neon.h>
#########################################################*/

/*
 gcc SMID flags
   - SSE is default gcc (no flag needed)
   - AVX is -mavx
   - AVX2 is -mavx2
*/

#ifndef VECTORWRAP_H
#define VECTORWRAP_H

/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOH: Start Of Header
'  o sec-01:
'    - Included libraries
'  o sec-02:
'    - variable defintions 
'  o sec-03:
'    - Structures
'  o sec-04:
'    - memBlock structure macros/functions (no functions)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-01:
^  - Included libraries
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifdef AVX512
   #include "vectWrapAVX512.h"

#elif defined AVX2
   #include "vectWrapAVX2.h"

#elif defined SSE4
   #define SSE2 /*Make sure not defined*/
   #include "vectWrapSSE.h"

#elif defined NEON
   #include "vectWrapNEON.h"

#elif defined NEON64
   #include "vectWrapNEON64.h"

/*SSE2 is last so it can be put as the fallback option*/
#elif defined SSE2
   #include "vectWrapSSE.h"

#else
   #define DEFNOVECT /*Not using vectors*/
   #define defVectBits 8
   #define defVectBytes 1
   #define defNum8BitElms 1
   #define defNum16BitElms 1
   #define defNum32BitElms 1

#endif

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-02:
^  - variable defintions 
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef DEFNOVECT /*If I am using vector support*/

#include <limits.h>
#include <stdlib.h>
/*stdint.h from vector libraries*/

/*Probably in stdint, but must makking sure I have the
`maximum value for each uintx_t value
*/
#define defMaxUI64 0xFFFFFFFFFFFFFFFF
#define defMaxUI32 0xFFFFFFFF
#define defMaxUI16 0xFFFF
#define defMaxUI8  0xFF

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-03:
^  - Structures
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

typedef struct vectMemBlock
{ /*vectMemBlock*/
   int8_t *alnPtr;  /*Holds aligned memory address*/
   int8_t *freePtr; /*Holds address to free*/
}vectMemBlock;

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec-04:
^  - memBlock structure macros/functions
^  o sec-04 macro-01: alnPointer
^    - returns an aligned pointer (as character pointer)
^  o sec-04 macro-02:
^    - Initialize a vectMemBlock structure (requires a
^      pointer to a vectMemBlock structuer)
^  o sec-04 macro-03:
^    - Make and initialize a memory block structure
^  o sec-04 macro-04:
^    - Allocate buffer with malloc for a vectMemBlock
^      structure makeVectMemBlock
^  o sec-04 macro-04:
^    - Allocate buffer with calloc for a vectMemBlock
^      structure makeVectMemBlock
^  o sec-04 macro-06:
^    - Free a vectMemBlock structure
^  o sec-04 macro-07:
^    - Free the buffer stored in a vectMemBlock structure
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*--------------------------------------------------------\
| Sec-04 Macro-01: alnPointer
|  - returns an aligned pointer (as character pointer)
\--------------------------------------------------------*/
#define alnPointer(alnPtr)(\
   (int8_t *)(\
        ((uintptr_t) (alnPtr))\
      & (~((uintptr_t) (defVectBytes - 1)))\
   )\
   /* Logic:
   `  defVectByes - 1
   `    gives 0xf (1111) for 16, 0xf1 (11111) for 32,
   `    and 0x4f (111111) for 64 
   ` ~ (defVectBytes - 1)
   `    gives 111...0000 (16 - 1), 111...00000 (32 - 1),
   `    or 111...000000 (64 - 1)
   ` freePtr & (~ (defVectByes - 1))
   `    Clears the first 4 (16), 5 (32), or 6 (64) bits
   `    in the pointer. Since eacy bit represents a byte
   `    this means that the pointer is now aligned for
   `    use with vectors.
   */\
)/*alnPointer*/

/*--------------------------------------------------------\
| Sec-04 Macro-02:
| Initialize a vectMemBlock structure (pionters are 0/null)
\--------------------------------------------------------*/
#define initVectMemBlock(memBlockST) \
{ \
   (memBlockST)->freePtr = 0;\
   (memBlockST)->alnPtr = 0;\
} /*initVectMemBlock*/

/*--------------------------------------------------------\
| Sec-04 Macro-03:
| Returns a initalized, allocated vectMemBlock structure
\--------------------------------------------------------*/
#define makeVectMemBlock() \
   (malloc(sizeof(struct vectMemBlock)))

/*--------------------------------------------------------\
| Sec-04 Macro-04:
| Allocate buffer with malloc for a vectMemBlock structure
\--------------------------------------------------------*/
#define mallocVectMemBlock(memBlockST, size)\
{ /*mallocVectMemBlock*/ \
   if((memBlockST)->freePtr != 0)\
   { /*If: there is already a buffer here*/\
      free((memBlockST)->freePtr);\
      initVectMemBlock((memBlockST));\
   } /*If: there is already a buffer here*/\
   \
   /*Allocate unaligned memory*/\
   (memBlockST)->freePtr = \
      malloc(sizeof(int8_t) * ((size) + defVectBytes));\
   \
   if((memBlockST)->freePtr != 0)\
      (memBlockST)->alnPtr=alnPointer(memBlockST->freePtr)\
} /*mallocVectMemBlock*/

/*--------------------------------------------------------\
| Sec-04 Macro-05:
| Allocate buffer with calloc for a vectMemBlock structure
\--------------------------------------------------------*/
#define callocVectMemBlock(memBlockST, size)\
{ /*callocVectMemBlock*/ \
   if((memBlockST)->freePtr != 0)\
   { /*If: there is already a buffer here*/\
      free((memBlockST)->freePtr);\
      initVectMemBlock((memBlockST));\
   } /*If: there is already a buffer here*/\
   \
   /*Allocate unaligned memory*/\
   (memBlockST)->freePtr = \
      calloc(((size) + defVectBytes), sizeof(int8_t));\
   \
   if((memBlockST)->freePtr != 0)\
      (memBlockST)->alnPtr=alnPointer(memBlockST->freePtr)\
} /*callocVectMemBlock*/

/*--------------------------------------------------------\
| Sec-04 Macro-06:
| Free a vectMemBlock structure
\--------------------------------------------------------*/
#define freeVectMemBlockBuff(memBlockST)\
{ /*freeMemBlockBuff*/ \
   if((memBlockST) != 0)\
   { /*If: I have a memory block to free a buffer from*/\
      if((memBlockST)->freePtr != 0)\
         free((memBlockST)->freePtr);\
      initVectMemBlock((memBlockST));\
   } /*If: I have a memory block to free a buffer from*/\
} /*freeMemBlockBuff*/ \

/*--------------------------------------------------------\
| Sec-04 Macro-07:
| Free the buffer stored in a vectMemBlock structure
\--------------------------------------------------------*/
#define freeVectMemBlock(memBlockST)\
{ /*freeMemBlock*/ \
   if((memBlockST) != 0)\
   { /*If: I have a memory block to free*/\
      freeVectMemBlockBuff(memBlockST)\
      free((memBlockST));\
      (memBlockST) = 0;\
   } /*If: I have a memory block to free*/\
} /*freeMemBlock*/ \

#endif /*If I am using vector support*/
#endif
