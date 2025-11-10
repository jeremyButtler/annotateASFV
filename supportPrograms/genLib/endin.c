/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' endin SOF: Start Of File
'   - converts 
'   o header:
'     - defined variables
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
|   - defined variables
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

/*these are assumptions that are not rigid, but if false
`  would be devisitating to the coding community. So,
`  though the number of bits in a byte is not formal it
`  is still very unlikely to be different on any modern
`  machine. Basically, it might as well be in the
`  standard.
*/
#define def_bitInByte_endin 8
#define def_ulHalf_endin (sizeof(unsigned long) << 2)
#define def_uiHalf_endin (sizeof(unsigned int) << 2)
#define def_usHalf_endin (sizeof(unsigned short) << 2)
   /*shifting by 2 instead of 3 to get half of the size
   `   of the datatype
   */

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
){
   #ifdef BIG_ENDIN
      return valUL; /*already in bigendin format*/
   #else
      unsigned char ucByte = 0;
      unsigned long retUL = 0;

      for(
         ucByte = def_bitInByte_endin;
         ucByte < ( sizeof(unsigned long) << 3 );
         ucByte += (def_bitInByte_endin << 1)
      ){ /*Loop: convert to bigendin*/
         
         retUL |=
            (
                (valUL << ucByte)
              & (
                      0xff
                   << (
                         def_ulHalf_endin
                       + (ucByte >> 1) - 4
                      )
                 )
            );
            /*Logic:
            `   - numBits: 32 + (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 32 is half of a unsigned long
            `     o 32 + ( 8 >> 1) - 4 = 32 +  4 - 4 = 32
            `     o 32 + (24 >> 1) - 4 = 32 + 12 - 4 = 40
            `     o 32 + (40 >> 1) - 4 = 32 + 20 - 4 = 48
            `     o 32 + (56 >> 1) - 4 = 32 + 28 - 4 = 56
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff00000000
            `     o 24: 0xff0000000000
            `     o 40: 0xff000000000000
            `     o 56: 0xff00000000000000
            `   - shift: valUL >> ucByte
            `     - shifts value to swap position
            `   - retUL |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */

         retUL |=
            (
                 (valUL >> ucByte)
               & (
                       0xff
                    << (
                           def_ulHalf_endin
                         - (ucByte >> 1) - 4
                       )
                  )
            );
            /*Logic:
            `   - numBits: 32 - (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 32 is half of a unsigned long
            `     o 32 - ( 8 >> 1) - 4 = 32 -  4 - 4 = 24
            `     o 32 - (24 >> 1) - 4 = 32 - 12 - 4 = 16
            `     o 32 - (40 >> 1) - 4 = 32 - 20 - 4 = 8
            `     o 32 - (56 >> 1) - 4 = 32 - 28 - 4 = 0
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff000000
            `     o 24: 0xff0000
            `     o 40: 0xff00
            `     o 56: 0xff
            `   - shift: valUL >> ucByte
            `     - shifts value to swap position
            `   - retUL |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */
      }  /*Loop: convert to bigendin*/

      return retUL;
   #endif
} /*ulToBig_endin*/

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
){
   #ifdef BIG_ENDIN
      return valUI; /*already in bigendin format*/
   #else
      unsigned char ucByte = 0;
      unsigned int retUI = 0;

      for(
         ucByte = def_bitInByte_endin;
         ucByte < ( sizeof(unsigned int) << 3 );
         ucByte += (def_bitInByte_endin << 1)
      ){ /*Loop: convert to bigendin*/
         
         retUI |=
            (
                (valUI << ucByte)
              & (
                      0xff
                   << (
                         def_uiHalf_endin
                       + (ucByte >> 1) - 4
                      )
                 )
            );
            /*Logic:
            `   - numBits: 16 + (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 16 is half of a unsigned int
            `     o 16 + ( 8 >> 1) - 4 = 16 +  4 - 4 = 16
            `     o 16 + (24 >> 1) - 4 = 16 + 12 - 4 = 24
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff0000
            `     o 24: 0xff000000
            `   - shift: valUI >> ucByte
            `     - shifts value to swap position
            `   - retUI |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */

         retUI |=
            (
                 (valUI >> ucByte)
               & (
                       0xff
                    << (
                           def_uiHalf_endin
                         - (ucByte >> 1) - 4
                       )
                  )
            );
            /*Logic:
            `   - numBits: 16 - (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 16 is half of a unsigned int
            `     o 16 - ( 8 >> 1) - 4 = 16 -  4 - 4 = 8
            `     o 16 - (24 >> 1) - 4 = 16 - 12 - 4 = 0
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff00
            `     o 24: 0xff
            `   - shift: valUI >> ucByte
            `     - shifts value to swap position
            `   - retUI |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */
      }  /*Loop: convert to bigendin*/

      return retUI;
   #endif
} /*uiToBig_endin*/

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
){
   #ifdef BIG_ENDIN
      return valUS; /*already in bigendin format*/
   #else
      unsigned char ucByte = 0;
      unsigned short retUS = 0;

      for(
         ucByte = def_bitInByte_endin;
         ucByte < ( sizeof(unsigned short) << 3 );
         ucByte += (def_bitInByte_endin << 1)
      ){ /*Loop: convert to bigendin*/
         
         retUS |=
            (
                (valUS << ucByte)
              & (
                      0xff
                   << (
                         def_usHalf_endin
                       + (ucByte >> 1) - 4
                      )
                 )
            );
            /*Logic:
            `   - numBits: 8 + (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 8 is half of a unsigned short
            `     o 8 + ( 8 >> 1) - 4 = 8 +  4 - 4 = 8
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff00
            `   - shift: valUS >> ucByte
            `     - shifts value to swap position
            `   - retUS |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */

         retUS |=
            (
                 (valUS >> ucByte)
               & (
                       0xff
                    << (
                           def_usHalf_endin
                         - (ucByte >> 1) - 4
                       )
                  )
            );
            /*Logic:
            `   - numBits: 8 - (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 8 is half of a unsigned short
            `     o 8 - ( 8 >> 1) - 4 = 8 -  4 - 4 = 0
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff
            `   - shift: valUS >> ucByte
            `     - shifts value to swap position
            `   - retUS |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */
      }  /*Loop: convert to bigendin*/

      return retUS;
   #endif
} /*usToBig_endin*/

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
){
   #ifdef BIG_ENDIN
      unsigned char ucByte = 0;
      unsigned long retUL = 0;

      for(
         ucByte = def_bitInByte_endin;
         ucByte < ( sizeof(unsigned long) << 3 );
         ucByte += (def_bitInByte_endin << 1)
      ){ /*Loop: convert to little endin*/
         
         retUL |=
            (
                (valUL << ucByte)
              & (
                      0xff
                   << (
                         def_ulHalf_endin
                       + (ucByte >> 1) - 4
                      )
                 )
            );
            /*Logic:
            `   - numBits: 32 + (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 32 is half of a unsigned long
            `     o 32 + ( 8 >> 1) - 4 = 32 +  4 - 4 = 32
            `     o 32 + (24 >> 1) - 4 = 32 + 12 - 4 = 40
            `     o 32 + (40 >> 1) - 4 = 32 + 20 - 4 = 48
            `     o 32 + (56 >> 1) - 4 = 32 + 28 - 4 = 56
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff00000000
            `     o 24: 0xff0000000000
            `     o 40: 0xff000000000000
            `     o 56: 0xff00000000000000
            `   - shift: valUL >> ucByte
            `     - shifts value to swap position
            `   - retUL |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */

         retUL |=
            (
                 (valUL >> ucByte)
               & (
                       0xff
                    << (
                           def_ulHalf_endin
                         - (ucByte >> 1) - 4
                       )
                  )
            );
            /*Logic:
            `   - numBits: 32 - (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 32 is half of a unsigned long
            `     o 32 - ( 8 >> 1) - 4 = 32 -  4 - 4 = 24
            `     o 32 - (24 >> 1) - 4 = 32 - 12 - 4 = 16
            `     o 32 - (40 >> 1) - 4 = 32 - 20 - 4 = 8
            `     o 32 - (56 >> 1) - 4 = 32 - 28 - 4 = 0
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff000000
            `     o 24: 0xff0000
            `     o 40: 0xff00
            `     o 56: 0xff
            `   - shift: valUL >> ucByte
            `     - shifts value to swap position
            `   - retUL |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */
      }  /*Loop: convert to little endin*/

      return retUL;
   #else
      return valUL; /*already in little endin format*/
   #endif
} /*ulToLittle_endin*/

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
){
   #ifdef BIG_ENDIN
      unsigned char ucByte = 0;
      unsigned int retUI = 0;

      for(
         ucByte = def_bitInByte_endin;
         ucByte < ( sizeof(unsigned int) << 3 );
         ucByte += (def_bitInByte_endin << 1)
      ){ /*Loop: convert to little endin*/
         
         retUI |=
            (
                (valUI << ucByte)
              & (
                      0xff
                   << (
                         def_uiHalf_endin
                       + (ucByte >> 1) - 4
                      )
                 )
            );
            /*Logic:
            `   - numBits: 16 + (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 16 is half of a unsigned int
            `     o 16 + ( 8 >> 1) - 4 = 16 +  4 - 4 = 16
            `     o 16 + (24 >> 1) - 4 = 16 + 12 - 4 = 24
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff0000
            `     o 24: 0xff000000
            `   - shift: valUI >> ucByte
            `     - shifts value to swap position
            `   - retUI |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */

         retUI |=
            (
                 (valUI >> ucByte)
               & (
                       0xff
                    << (
                           def_uiHalf_endin
                         - (ucByte >> 1) - 4
                       )
                  )
            );
            /*Logic:
            `   - numBits: 16 - (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 16 is half of a unsigned int
            `     o 16 - ( 8 >> 1) - 4 = 16 -  4 - 4 = 8
            `     o 16 - (24 >> 1) - 4 = 16 - 12 - 4 = 0
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff00
            `     o 24: 0xff
            `   - shift: valUI >> ucByte
            `     - shifts value to swap position
            `   - retUI |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */
      }  /*Loop: convert to little endin*/

      return retUI;
   #else
      return valUI; /*already in little endin format*/
   #endif
} /*uiToLittle_endin*/

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
){
   #ifdef BIG_ENDIN
      unsigned char ucByte = 0;
      unsigned short retUS = 0;

      for(
         ucByte = def_bitInByte_endin;
         ucByte < ( sizeof(unsigned short) << 3 );
         ucByte += (def_bitInByte_endin << 1)
      ){ /*Loop: convert to little endin*/
         
         retUS |=
            (
                (valUS << ucByte)
              & (
                      0xff
                   << (
                         def_usHalf_endin
                       + (ucByte >> 1) - 4
                      )
                 )
            );
            /*Logic:
            `   - numBits: 8 + (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 8 is half of a unsigned short
            `     o 8 + ( 8 >> 1) - 4 = 8 +  4 - 4 = 8
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff00
            `   - shift: valUS >> ucByte
            `     - shifts value to swap position
            `   - retUS |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */

         retUS |=
            (
                 (valUS >> ucByte)
               & (
                       0xff
                    << (
                           def_usHalf_endin
                         - (ucByte >> 1) - 4
                       )
                  )
            );
            /*Logic:
            `   - numBits: 8 - (ucByte >> 1) - 4
            `     - number of bits need to shift 0xff with
            `       to clear everything except the target
            `       - 8 is half of a unsigned short
            `     o 8 - ( 8 >> 1) - 4 = 8 -  4 - 4 = 0
            `   - mask: 0xff << numBits
            `     - mask to clear everything except swapped
            `       value
            `     o 8:  0xff
            `   - shift: valUS >> ucByte
            `     - shifts value to swap position
            `   - retUS |= (shift & mask)
            `     - clears bits not in swap & adds to
            `       return value
            */
      }  /*Loop: convert to little endin*/

      return retUS;
   #else
      return valUS; /*already in little endin format*/
   #endif
} /*usToLittle_endin*/

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
){
   unsigned char ucByte = 0;
   unsigned long retUL = 0;

   for(
      ucByte = def_bitInByte_endin;
      ucByte < ( sizeof(unsigned long) << 3 );
      ucByte += (def_bitInByte_endin << 1)
   ){ /*Loop: change endin*/
      
      retUL |=
         (
             (valUL << ucByte)
           & (
                   0xff
                << (
                      def_ulHalf_endin
                    + (ucByte >> 1) - 4
                   )
              )
         );
         /*Logic:
         `   - numBits: 32 + (ucByte >> 1) - 4
         `     - number of bits need to shift 0xff with
         `       to clear everything except the target
         `       - 32 is half of a unsigned long
         `     o 32 + ( 8 >> 1) - 4 = 32 +  4 - 4 = 32
         `     o 32 + (24 >> 1) - 4 = 32 + 12 - 4 = 40
         `     o 32 + (40 >> 1) - 4 = 32 + 20 - 4 = 48
         `     o 32 + (56 >> 1) - 4 = 32 + 28 - 4 = 56
         `   - mask: 0xff << numBits
         `     - mask to clear everything except swapped
         `       value
         `     o 8:  0xff00000000
         `     o 24: 0xff0000000000
         `     o 40: 0xff000000000000
         `     o 56: 0xff00000000000000
         `   - shift: valUL >> ucByte
         `     - shifts value to swap position
         `   - retUL |= (shift & mask)
         `     - clears bits not in swap & adds to
         `       return value
         */

      retUL |=
         (
              (valUL >> ucByte)
            & (
                    0xff
                 << (
                        def_ulHalf_endin
                      - (ucByte >> 1) - 4
                    )
               )
         );
         /*Logic:
         `   - numBits: 32 - (ucByte >> 1) - 4
         `     - number of bits need to shift 0xff with
         `       to clear everything except the target
         `       - 32 is half of a unsigned long
         `     o 32 - ( 8 >> 1) - 4 = 32 -  4 - 4 = 24
         `     o 32 - (24 >> 1) - 4 = 32 - 12 - 4 = 16
         `     o 32 - (40 >> 1) - 4 = 32 - 20 - 4 = 8
         `     o 32 - (56 >> 1) - 4 = 32 - 28 - 4 = 0
         `   - mask: 0xff << numBits
         `     - mask to clear everything except swapped
         `       value
         `     o 8:  0xff000000
         `     o 24: 0xff0000
         `     o 40: 0xff00
         `     o 56: 0xff
         `   - shift: valUL >> ucByte
         `     - shifts value to swap position
         `   - retUL |= (shift & mask)
         `     - clears bits not in swap & adds to
         `       return value
         */
   }  /*Loop: change endin*/

   return retUL;
} /*ulSwap_endin*/

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
){
   unsigned char ucByte = 0;
   unsigned int retUI = 0;

   for(
      ucByte = def_bitInByte_endin;
      ucByte < ( sizeof(unsigned int) << 3 );
      ucByte += (def_bitInByte_endin << 1)
   ){ /*Loop: swap to oppisite endin*/
      
      retUI |=
         (
             (valUI << ucByte)
           & (
                   0xff
                << (
                      def_uiHalf_endin
                    + (ucByte >> 1) - 4
                   )
              )
         );
         /*Logic:
         `   - numBits: 16 + (ucByte >> 1) - 4
         `     - number of bits need to shift 0xff with
         `       to clear everything except the target
         `       - 16 is half of a unsigned int
         `     o 16 + ( 8 >> 1) - 4 = 16 +  4 - 4 = 16
         `     o 16 + (24 >> 1) - 4 = 16 + 12 - 4 = 24
         `   - mask: 0xff << numBits
         `     - mask to clear everything except swapped
         `       value
         `     o 8:  0xff0000
         `     o 24: 0xff000000
         `   - shift: valUI >> ucByte
         `     - shifts value to swap position
         `   - retUI |= (shift & mask)
         `     - clears bits not in swap & adds to
         `       return value
         */

      retUI |=
         (
              (valUI >> ucByte)
            & (
                    0xff
                 << (
                        def_uiHalf_endin
                      - (ucByte >> 1) - 4
                    )
               )
         );
         /*Logic:
         `   - numBits: 16 - (ucByte >> 1) - 4
         `     - number of bits need to shift 0xff with
         `       to clear everything except the target
         `       - 16 is half of a unsigned int
         `     o 16 - ( 8 >> 1) - 4 = 16 -  4 - 4 = 8
         `     o 16 - (24 >> 1) - 4 = 16 - 12 - 4 = 0
         `   - mask: 0xff << numBits
         `     - mask to clear everything except swapped
         `       value
         `     o 8:  0xff00
         `     o 24: 0xff
         `   - shift: valUI >> ucByte
         `     - shifts value to swap position
         `   - retUI |= (shift & mask)
         `     - clears bits not in swap & adds to
         `       return value
         */
   }  /*Loop: swap to oppisite endin*/

   return retUI;
} /*uiSwap_endin*/

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
){
   unsigned char ucByte = 0;
   unsigned short retUS = 0;

   for(
      ucByte = def_bitInByte_endin;
      ucByte < ( sizeof(unsigned short) << 3 );
      ucByte += (def_bitInByte_endin << 1)
   ){ /*Loop: swap to oppisite endin*/
      
      retUS |=
         (
             (valUS << ucByte)
           & (
                   0xff
                << (
                      def_usHalf_endin
                    + (ucByte >> 1) - 4
                   )
              )
         );
         /*Logic:
         `   - numBits: 8 + (ucByte >> 1) - 4
         `     - number of bits need to shift 0xff with
         `       to clear everything except the target
         `       - 8 is half of a unsigned short
         `     o 8 + ( 8 >> 1) - 4 = 8 +  4 - 4 = 8
         `   - mask: 0xff << numBits
         `     - mask to clear everything except swapped
         `       value
         `     o 8:  0xff00
         `   - shift: valUS >> ucByte
         `     - shifts value to swap position
         `   - retUS |= (shift & mask)
         `     - clears bits not in swap & adds to
         `       return value
         */

      retUS |=
         (
              (valUS >> ucByte)
            & (
                    0xff
                 << (
                        def_usHalf_endin
                      - (ucByte >> 1) - 4
                    )
               )
         );
         /*Logic:
         `   - numBits: 8 - (ucByte >> 1) - 4
         `     - number of bits need to shift 0xff with
         `       to clear everything except the target
         `       - 8 is half of a unsigned short
         `     o 8 - ( 8 >> 1) - 4 = 8 -  4 - 4 = 0
         `   - mask: 0xff << numBits
         `     - mask to clear everything except swapped
         `       value
         `     o 8:  0xff
         `   - shift: valUS >> ucByte
         `     - shifts value to swap position
         `   - retUS |= (shift & mask)
         `     - clears bits not in swap & adds to
         `       return value
         */
   }  /*Loop: swap to oppisite endin*/

   return retUS;
} /*usSwap_endin*/

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
