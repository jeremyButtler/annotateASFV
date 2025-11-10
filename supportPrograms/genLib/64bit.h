/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' 64bit SOF:
'   - holds data type and macro to get the maximum long
'     on a system
'   - this will not always be a 64 bit, but the values
'     will still behave the same so long as bytes are 8
'     bits
'   o .h var01:
'     - define ul_64bit and sl_64bit variable
'   o .h def01:
'     - definitions dealing with size of the ulong type
'   o .h def02:
'     - values to change byte lengths to get ulong lengths
'   o .h fun01: emptyBytes_64bit
'     - find number of empty bytes in a ulong
'   o .h fun02: bytesInUL_64bit
'     - find number of bytes in a single ulong
'   o .h fun03: byteLenToULLen_64bit
'     - converts length in bytes to number of ulongs
'       needed to hold the bytes
'   o .h fun04: byteLenToFullULLen_64bit
'      - converts length in bytes to number of fully
'        filled ulongs (a paritially filled ulong is lost)
'   o .h fun05: rshiftByte_64bit
'     - shift a byte right by x bytes; over shifts go to 0
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Var01:
|   - define ul_64bit and sl_64bit variable
\-------------------------------------------------------*/
#ifdef NOUL
   typedef unsigned char ul_64bit;
   typedef signed char sl_64bit;
#else
   #ifdef PLAN9_64
      /*plan9 defines long as 32 bit, so need long long*/
      typedef unsigned long long ul_64bit;
      typedef signed long long sl_64bit;
   #else
      typedef unsigned long ul_64bit;
      typedef signed long sl_64bit;
   #endif
#endif

/*-------------------------------------------------------\
| Def01:
|   - definitions dealing with size of the ulong type
\-------------------------------------------------------*/
#define def_bitsPerChar_64bit 8 /*number bits / byte*/
#define def_shiftBitsToBytes_64bit 3
   /*do `variant >> def_shiftBitsToBytes_64bit` to get
   `  the number of bytes from bits
   `or do `variant << def_shiftBitsToBytes_64bit` to get
   `  the number of bits from bytes
   */

#define def_bytesInUL_64bit (sizeof(ul_64bit))
   /*number of bytes in the long*/
#define def_bitsInUL_64bit (def_bytesInUL_64bit << def_shiftBitsToBytes_64bit)
   /*number of bits in a long*/

#define def_ulNegBit_64bit ( ((ul_64bit) -1) & (((ul_64bit) -1) >> 1) )
   /*negative bit in a long*/
#define def_maxUL_64bit ((ul_64bit) -1)
   /*maximum value in ulong*/
#define def_maxSL_64bit ( (sl_64bit) ((ul_64bit) ~def_ulNegBit_64bit) )
   /*maximum value in ulong*/
#define def_minSL_64bit ((sl_64bit) def_ulNegBit_64bit)
   /*maximum value in ulong*/

/*-------------------------------------------------------\
| Def02:
|   - values to change byte lengths to get ulong lengths
\-------------------------------------------------------*/

/*do `byte_length & def_modUL_64bit` to find the number of
`   trailing bytes in the last ulong
` do `byte_length - (byte_length & def_modUL_64bit)` to
`   find the number of empty bytes in a ulong
*/
#define def_modUL_64bit (sizeof(ul_64bit ) - 1)
   /*The mod values are 7 (longs), 3 (ints), and 1
   `  (shorts), so, 1 - number of bytes will give me the
   `  mod value
   */

/*Used to convert a length in bytes to length in ulongs
`This should be compliled in.
`  This will only work for up to 128 bit numbers, after
`  that it will fail. However, it is very unlikely a
`  long will be even 128 bits.
*/
#define def_byteLenToULLen_64bit ( ((sizeof(ul_64bit) >> 1) - (sizeof(ul_64bit) >> 3)) - ((sizeof(ul_64bit) >> 3) * (sizeof(ul_64bit) >> 4)) )

/*Logic:
`   - 32bitSize: sizeof(ul_64bit) >> 1:
`     o This gets the number of bytes to shift to get
`       an 32 bit or 16 bit number to one bytes. Basically
`       it is divide by 2
`   - (sizeof(ul_64bit) >> 3)
`     o This adds in 1 if I have eight bytes. It is needed
`       as the 64 bit correction step
`   - 64bitSize: 32bitSize - (sizeof(ul_64bit) >> 3
`     - (sizeof(ul_64bit) >> 1) - (sizeof(ul_64bit) >> 3)
`     o Gives me the number of bits I need to shift to
`       get the number of characters in an long
`   - 128bitSize: 64bitSize - (sizeof(ul_64bit) >> 3) * (sizeof(ul_64bit) >> 4))
`     - (sizeof(ul_64bit >> 1) - sizeof(ul_64bit) >> 3) - ((sizeof(ul_64bit) >> 3) * (sizeof(ul_64bit) >> 4))
`     o used to correct for a 128 bit number.
`   - for values > 128 bit you need to provided a new
`     correction step
*/

/*-------------------------------------------------------\
| Fun01: emptyBytes_64bit
|   - find number of empty bytes in a ulong
| Input:
|   - lenMacroVar1:
|     o length of string in bytes
| Output:
|   - Returns:
|     o number of extra bytes at the end
\-------------------------------------------------------*/
#define emptyBytes_64bit(lenMacroVar1) ((lenMacroVar1) - (def_modUL_64bit & (lenMacroVar1)) )

/*-------------------------------------------------------\
| Fun02: bytesInUL_64bit
|   - find number of bytes in a single ulong
| Input:
|   - lenMacroVar2:
|     o length of string in bytes
| Output:
|   - Returns:
|     o number of extra bytes at the end
\-------------------------------------------------------*/
#define bytesInUL_64bit(lenMacroVar2) (def_bytesInUL_64bit - emptyBytes_64bit((lenMacroVar2)))

/*-------------------------------------------------------\
| Fun03: byteLenToULLen_64bit
|   - converts length in bytes to number of ulongs needed
|     to hold the bytes
| Input:
|   - lenMacroVar3:
|     o number of bytes (byte length)
| Output:
|   - Returns:
|     o number of ulongs needed to hold lenSI bytes
\-------------------------------------------------------*/
#define byteLenToULLen_64bit(lenMacroVar3) ( ((lenMacroVar3) >> def_byteLenToULLen_64bit) + (emptyBytes_64bit((lenMacroVar3)) > 0) )

/*-------------------------------------------------------\
| Fun04: byteLenToFullULLen_64bit
|   - converts length in bytes to number of fully filled
|     ulongs (a paritially filled ulong is lost)
| Input:
|   - lenMacroVar4:
|     o number of bytes (byte length)
| Output:
|   - Returns:
|     o number of ulongs needed to hold lenSI bytes
\-------------------------------------------------------*/
#define byteLenToFullULLen_64bit(lenMacroVar4) ((lenMacroVar4) >> def_byteLenToULLen_64bit)

/*-------------------------------------------------------\
| Fun05: rshiftByte_64bit
|   - shift a byte right by x bytes (over shifts go to 0)
| Input:
|   - val:
|     o value to shift right
|   - byte:
|     o number of bytes to shift by
| Output:
|   - Returns:
|     o shifted value
|     o 0 if over shifted
\-------------------------------------------------------*/
#define rshiftByte_64bit(val, byte) ( ( ((ul_64bit) (val)) << ( ( ((ul_64bit) (byte)) & (-(( ((signed int) sizeof(ul_64bit)) - ((signed int) (byte)) ) > 0 ))) << 3 ) ) & (-(( ((signed int) sizeof(ul_64bit)) - ((signed int) (byte)) ) > 0 )))
/* Logic:
`    - right shift by bytes
`    - boundsCheck: (sizeof(ul_64bit)) - byte) ) > 0:
`      o check if shift is out of bounds
`      o 0 if out of bounds (ul_64bit - byte) <= 0
`      o 1 if in bounds of ul_64bit
`    - boundsMask: - boundsCheck: 
`      o 0 if out of bounds
`      o -1 if in of bounds
`    - getBytes: byte & boundsMask:
`      o check if keep bytes
`      o 0 if out of bounds
`      o bytes if in of bounds
`    - shiftVal: getBytes << 3:
`      o get the number of bits to shift
`      o 0 if out of bounds
`      0 number bytes * 8
`    - shift: val << shiftVal:
`      o shift value to correct position
`      0 value if out of bounds
`      0 shifted value
`    - finalVal: shift & boundsMask:
`      o sets out of bounds values to 0
`      o 0 if out of bounds
`      o shifted value if in bounds
*/
