/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' strFun: SOF (Start Of File)
'   - general non-optimized string functions
'   o header:
'     - included libraries
'   o fun01: pad_strFun
'     - pads an input c-string to the desired length
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
#endif

#include "strFun.h"

/*-------------------------------------------------------\
| Fun01: pad_strFun
|   - pads an input c-string to the desired length
| Input:
|   - inStr:
|     o c-string to add padding to (should have room)
|   - padSC:
|     o character to pad with
|   - lenSI:
|     o length string should be when done padding
|   - dirSC:
|     o direction to apply padding in
|       * 0: left padding
|       * 1: center padding
|       * 2: right padding
| Output:
|   - Modifies:
|     o inStr to have padding
\-------------------------------------------------------*/
void
pad_strFun(
   signed char *inStr, /*c-string to pad*/
   signed char padSC,  /*character to pad with*/
   signed int lenSI,   /*final length after padding*/
   signed char dirSC,  /*direction to pad*/
){
   signed char *dupStr = inStr;
   signed char *cpStr = inStr;

   signed int padLenSI = 0;
   signed long inLenSL = 0;

   while(inStr[inLenSL])
      ++inLenSL;

   cpStr += inLenSL;

   if(inLenSL > lenSI)
      return;
   else
      padLenSI = lenSI - inLenSL;

   if(dirSC == 1)
      padLenSI >>= 1;


   if(dirSC == 2)
   { /*Else If: right padding*/
      while(padLenSI--)
         inStr[inLenSL++] = padSC;
   } /*Else If: right padding*/

   else
   { /*Else: left or center padding*/
      dupStr = inLenSL + padLenSI;
      dupStr[1] = 0;

      while(*cpStr >= inStr)
         *dupStr-- = *cpStr++;

      cpStr = inStr + inLenSL + padLenSI;
      dupStr = inStr;
      lenSI -= padLenSI;

      while(padLenSI--)
         *dupStr++ = padSC;

      if(padSC == 1)
      { /*If: center padding*/
          while(lenSI--)
             *cpStr++ = padSC;
          *cpStr = 0;
      } /*If: center padding*/
   } /*Else: left or center padding*/
} /*pad_strFun*/
