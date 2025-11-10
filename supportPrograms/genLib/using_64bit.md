# Goal:

Describe the 64 bit .h library and how to use it ot
  leverage longs on a system.

File: 64bit.h

# 64bit overiew

It always assigns a 64 bit value to its long on a OS/CPU
  that is at least 64 bit. For non-64 bit OS's/CPU's it
  assigns the largest datatype possible.

The data types are `ul_64bit` for unsigned long
  and `sl_64bit` for signed long.

The functions and defined values should updated based on
  the system they are compiled on.

## using in your .h file

Not recomended.

The plan9 native compilers enforce no `#include` in .h
  files. So, the best way I found to use the 64bit.h
  library in .h files is to include this is your .h file
  header is to forward declare a new variable of the same
  size. The use that variable in your function
  declarations.

```
#ifdef NOUL
   typedef unsigned char ul_<file_name>
   typedef signed char sl_<file_name>
#else
   #ifdef PLAN9_64
      typedef unsigned long long ul_<file_name>
      typedef signed long long sl_<file_name>
   #else
      typedef unsigned long ul_<file_name>
      typedef signed long sl_<file_name>
   #endif
#endif
```

Relace the `<file_name>` with your own unique name. I
  would suggest your file name so it stays unique between
  files.

## Functions

### Converting byte lengths to long lengths

You can use byteLenToFullULLen_64bit to convert a length
  to the number of full longs it can have. For example,
  80 would go to ten longs. However, this will trim off
  the incomplete  long at the end, so ten would got to one
  long (not two).

- Example: `int newLen = byteLenToFullULLen_64bit(80)`
  - gives 10
- Example: `int newLen = byteLenToFullULLen_64bit(10)`
  - gives 1

You can use the byteLenToULLen_64bit to get the number of
  longs needed to store a length of bytes. This function
  will also return the long at the end.

- Example: `int newLen = byteLenToFullULLen_64bit(80)`
  - gives 10
- Example: `int newLen = byteLenToFullULLen_64bit(10)`
  - gives 2

### Finding number of bytes in longs

You can get the number of bytes in a single long, or in
  the last long (if length is more then one long) with the
  `bytesInUL_64bit()` function. You can get the inverse
  (number of bytes open) in the last long with
  the `emptyBytes_64bit()` function.

For both fuctions the input is the number of bytes in the
  long or long array.

- Example `int bytesInLastLong = bytesInUL_64bit(80)`
  - will give 8 (max possible) value (last long is full)
- Example `int bytesInLastLong = bytesInUL_64bit(79)`
  - will give 7 value (last long has one byte left)
- Example `int bytesInLastLong = bytesInUL_64bit(78)`
  - will give 6 value (last long has 2 bytes left)
- Example `int bytesInLastLong = bytesInUL_64bit(77)`
  - will give 5 value (last long has 3 bytes left)
- Example `int bytesInLastLong = bytesInUL_64bit(76)`
  - will give 4 value (last long has 4 bytes left)
- Example `int bytesInLastLong = bytesInUL_64bit(75)`
  - will give 3 value (last long has 5 bytes left)
- Example `int bytesInLastLong = bytesInUL_64bit(74)`
  - will give 2 value (last long has 6 bytes left)
- Example `int bytesInLastLong = bytesInUL_64bit(73)`
  - will give 1 value (last long has 7 bytes left)
- Example `int bytesInLastLong = bytesInUL_64bit(72)`
  - will give 8 (max possible) value (last long is full)

## Predefiend values

- hardcoded values, but should never change
  - `def_bitsPerChar_64bit`: number bits in one byte
    - hardcoded, so does not update
  - `def_shiftBitsToBytes_64bita: value to right shift bit
    count (`>>`) by to get bytes.
    - can left shift (`<<`) bytes to get bits
    - hardcoded, so does not update
- size of `ul_64bit` and `sl_64bit` data typesj
  - `def_bytesInUL_64bit`: number of bytes in a
    `ul_64bit` or a `sl_64bit` data type
  - `def_bitsInUL_64bit`: number of bits in `ul_64bit`
    (for slong, do `def_bitsInUL_64bit - 1`)
- minimum and maximum values
  - `def_maxUL_64bit`: maximum value for an unsigned long
  - `def_maxSL_64bit`: maximum value for a signed long
  - `def_minSL_64bit`: minimum value for a signed long
  - `def_ulNegBit_64bit`: position of the negative bit in
    a signed long
- length values
  - `def_modUL_64bit`: allows you to find how many bytes
    are in the last long from the length
    - `length_in_bytes & def_modUL_64bit` will give 1 to
      7 for not full and 0 for full
  - `def_byteLenToULLen_64bit` right shift by this to get
    the number of full longs (partially filled long at end
    is ignored)
    - `byte_length >> def_byteLenToULLen_64bit`
