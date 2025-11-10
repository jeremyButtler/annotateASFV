/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' strFun: SOF (Start Of File)
'   - general non-optimized string functions
'   o header:
'     - guards
'   o fun01: pad_strFun
'     - pads an input c-string to the desired length
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef GENERAL_STRING_FUNCTIONS_H
#define GENERAL_STRING_FUNCTIONS_H

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
);

#endif
