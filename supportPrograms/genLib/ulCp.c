/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' ulCp SOF: Start Of File
'   - uses longs to copy contents of a string
'   - these functions are slower than strcpy and strlen,
'     but they allow deliminators to be used.
'   - DOES NOT WORK ON **BIG ENDIN** SYSTEM
'   o header:
'     - included libraries
'   deliminator detection/making:
'     o .h fun01: rshiftByte_ulCp
'       - shift byte right by x bytes; over shifts go to 0
'     o .h fun02: mkDelim_ulCp
'       - make ulong_ulCp delimintor from character
'     o .h fun03: ifDelim_ulCp
'       - check if have deliminator in long
'     o .h fun04: ifNullDelim_ulCp
'       - checks for null or deliminator in ulong_ulCp
'     o .h fun05: ifNull_ulCp
'       - checks to see if null in ulong_ulCp
'     o .h fun06: ifLineUnix_ulCp
'       - checks if unix line end ('\n' or '\0')
'     o .h fun07: ifEndLine_ulCp
'       - checks if input long is end of line
'     o .h fun08: ifWhite_ulCp
'       - check if input long has white space
'         (' ', '\r', '\n', '\t', '\0')
'   String copying:
'     o fun09: cpLen_ulCp
'       - copies cpStr into dupStr using ulong_ulCps
'     o fun10: cpDelim_ulCp
'       - copies string until deliminator is found
'     o fun11: cpStr_ulCp
'       - copies string until \0
'     o fun12: cpLineUnix_ulCp
'       - copies string until end of line (\0, \n)
'     o fun13: cpLine_ulCp
'       - copies string until end of line (\0, \r, \n)
'       - ingores unused ascii characters (> 32, not '\t')
'     o fun14: cpWhite_ulCp
'       - copies string until white space
'   String length:
'     o fun15: lenStr_ulCp
'       - finds the length of a string using ulong_ulCps
'     o fun16: lenStrNull_ulCp
'       - finds the length of a string using ulong_ulCps
'         but also stops at null '\0'
'     o fun17: endStr_ulCp
'       - finds the end of a c-string ('\0')
'     o fun18: endLineUnix_ulCp
'      - finds the end of a c-string. This assumes that
'        the line ends in '\0' or '\n'
'     o fun19: endLine_ulCp
'       - finds the end of a c-string (all OS's; looks for
'         '\0', '\n', and '\r')
'     o fun20: endWhite_ulCp
'       - finds 1st white space ('\t',' ','\n','\0','\r')
'   String comparision:
'     o fun21: eql_ulCp
'       - compares two strings until deliminator is found
'     o fun22: eqlNull_ulCp
'       - compares two strings until null is found
'     o fun23: eqlWhite_ulCp
'       - compares two strings until white space is found
'   String cleanup:
'     o fun24: rmWhite_ulCp
'       - removes white space from c-string
'   String swap:
'      o fun25: swapDelim_ulCp
'        - swaps two strings until deliminator is found
'      o fun26: swapNull_ulCp
'        - swaps two strings until null
'      o fun27: shift_ulCp
'        - shifts a substring in a string up or down
'        - not tested, likely does not work
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
#endif

#include "ulCp.h"

/*-------------------------------------------------------\
| Fun09: cpLen_ulCp
|   - copies cpStr into dupStr using ulong_ulCps
| Input:
|   - dupStr:
|     o pointer to string to copy cpStr into
|   - cpStr:
|     o pointer to string to copy
|   - lenUI:
|     o length of cpStr
| Output:
|   - Modifies:
|     o  dupStr to hold lenUI characters from cpStr
\-------------------------------------------------------*/
#ifdef NOUL
   void
   cpLen_ulCp(
      signed char *dupStr,
      signed char *cpStr,
      unsigned int lenUI
   ){
      while(lenUI--)
         *dupStr++ = *cpStr++;
      *dupStr = '\0';
   } /*cpLen_ulCp*/
#else
   void
   cpLen_ulCp(
      signed char *dupStr,
      signed char *cpStr,
      unsigned int lenUI
   ){
      ulong_ulCp *cpUL = (ulong_ulCp *) cpStr;
      ulong_ulCp *dupUL = (ulong_ulCp *) dupStr;
      unsigned int uiChar = 0;
   
      for(
         uiChar= 0;
         uiChar < (lenUI >> def_shiftULBy_ulCp);
         ++uiChar
      ) dupUL[uiChar] = cpUL[uiChar];
   
      for(
         uiChar = ( lenUI - (lenUI & def_modUL_ulCp) );
         uiChar < lenUI;
        ++uiChar
      ) dupStr[uiChar] = cpStr[uiChar];
   
      dupStr[lenUI] = '\0';
   } /*cpLen_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun10: cpDelim_ulCp
|   - copies string until deliminator is found
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|     o this will give an unused variable warning on
|       -DNOUL, but allows user to turn off ulCp
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   cpDelim_ulCp(
      signed char *dupStr,
      signed char *cpStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      signed char *endStr = dupStr;
   
      while(*cpStr != delimSC)
         *endStr++ = *cpStr++;
      *endStr = '\0';
      return endStr - dupStr;
   } /*cpDelim_ulCp*/
#else
   unsigned int
   cpDelim_ulCp(
      signed char *dupStr,
      signed char *cpStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
      ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
      signed char *dupTmpStr = 0;
   
      while( ! ifDelim_ulCp(*cpUL, delimUL) )
         *dupUL++ = *cpUL++;
   
      cpStr = (signed char *) cpUL;
      dupTmpStr = (signed char *) dupUL;
   
      while( *cpStr != delimSC )
         *dupTmpStr++ = *cpStr++;
   
      *dupTmpStr = '\0';
      return dupTmpStr - dupStr;
   } /*cpDelim_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun11: cpStr_ulCp
|   - copies string until \0
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   cpStr_ulCp(
      signed char *dupStr,
      signed char *cpStr
   ){
      signed char *endStr = dupStr;
   
      while(*cpStr)
         *endStr++ = *cpStr++;
      *endStr = '\0';
      return endStr - dupStr;
   } /*cpStr_ulCp*/
#else
   unsigned int
   cpStr_ulCp(
      signed char *dupStr,
      signed char *cpStr
   ){
      ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
      ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
      signed char *dupTmpStr = 0;
   
      while( ! ifNull_ulCp(*cpUL) )
         *dupUL++ = *cpUL++;
   
      cpStr = (signed char *) cpUL;
      dupTmpStr = (signed char *) dupUL;
   
      while( *cpStr != '\0' )
         *dupTmpStr++ = *cpStr++;
   
      *dupTmpStr = '\0';
      return dupTmpStr - dupStr; /*number of char copied*/
   } /*cpStr_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun12: cpLineUnix_ulCp
|   - copies string until end of line (\0, \n)
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   cpLineUnix_ulCp(
      signed char *dupStr,
      signed char *cpStr
   ){
      signed char *endStr = dupStr;
   
      while( *cpStr & (~ '\n') )
         *endStr++ = *cpStr++;
      *endStr = '\0';
      return endStr - dupStr;
   } /*cpStr_ulCp*/
#else
   unsigned int
   cpLineUnix_ulCp(
      signed char *dupStr,
      signed char *cpStr
   ){
      ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
      ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
      signed char *dupTmpStr = 0;
   
      while( ! ifLineUnix_ulCp(*cpUL) )
         *dupUL++ = *cpUL++;
   
      cpStr = (signed char *) cpUL;
      dupTmpStr = (signed char *) dupUL;
   
      while( *cpStr & (~ '\n') )
         *dupTmpStr++ = *cpStr++;
   
      *dupTmpStr = '\0';
      return dupTmpStr - dupStr; /*number of char copied*/
   } /*cpLineUnix_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun13: cpLine_ulCp
|   - copies string until end of line (\0, \r, \n)
|   - ingores all unused ascii characters (> 32, not '\t')
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   cpLine_ulCp(
      signed char *dupStr,
      signed char *cpStr
   ){
      signed char *endStr = dupStr;
   
      while(
            *cpStr > '\r'  /*\r > \t > \n > \0*/
         || *cpStr == '\t' /*so catch tab case*/
      ) *endStr++ = *cpStr++;

      *endStr = '\0';
      return endStr - dupStr;
   } /*cpLine_ulCp*/
#else
   unsigned int
   cpLine_ulCp(
      signed char *dupStr,
      signed char *cpStr
   ){
      ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
      ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
      signed char *dupTmpStr = 0;
   
      while( ! ifEndLine_ulCp(*cpUL) )
         *dupUL++ = *cpUL++;
   
      cpStr = (signed char *) cpUL;
      dupTmpStr = (signed char *) dupUL;
   
      while(
            *cpStr > '\r'  /*\r > \t > \n > \0*/
         || *cpStr == '\t' /*so catch tab case*/
      ) *dupTmpStr++ = *cpStr++;
   
      *dupTmpStr = '\0';
      return dupTmpStr - dupStr; /*number of char copied*/
   } /*cpLine_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun14: cpWhite_ulCp
|   - copies string until white space
| Input:
|   - dupStr:
|     o Pointer to string to copy cpStr into
|   - cpStr:
|     o Pointer to string to copy
| Output:
|   - Modifies:
|     o  dupStr to hold the characters from cpStr
| Note:
|   - This will likely not be very good at copying short
|     strings.
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   cpWhite_ulCp(
      signed char *dupStr,
      signed char *cpStr
   ){

      signed char *endStr = dupStr;
   
      while(*cpStr > 32)
         *endStr++ = *cpStr++;
      *endStr = '\0';
      return endStr - dupStr;
   } /*cpWhite_ulCp*/
#else
   unsigned int
   cpWhite_ulCp(
      signed char *dupStr,
      signed char *cpStr
   ){
      ulong_ulCp *cpUL = (ulong_ulCp *) (cpStr);
      ulong_ulCp *dupUL = (ulong_ulCp *) (dupStr);
      signed char *dupTmpStr = 0;
   
      while( ! ifWhite_ulCp(*cpUL) )
         *dupUL++ = *cpUL++;
   
      cpStr = (signed char *) cpUL;
      dupTmpStr = (signed char *) dupUL;
   
      while( *cpStr > 33 )
         *dupTmpStr++ = *cpStr++;
   
      *dupTmpStr = '\0';
      return dupTmpStr - dupStr; /*number of char copied*/
   } /*cpWhite_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun15: lenStr_ulCp
|   - finds the length of a string using ulong_ulCps
| Input:
|   - inStr:
|     o c-string or string with deliminator to find length
|       of
|   - delimUL:
|     o deliminator (as ulong_ulCp (fun02)) at end of
|       string
|     o this will give an unused variable warning on
|       -DNOUL, but allows user to turn off ulCp
|   - delimSC:
|     o deliminator (as char) at end of string
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   lenStr_ulCp(
      signed char *inStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      signed char *tmpStr = inStr;
   
      while(*tmpStr != delimSC)
         ++tmpStr;
   
      return tmpStr - inStr;
   } /*lenStr_ulCp*/
#else
   unsigned int
   lenStr_ulCp(
      signed char *inStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
      signed char *tmpStr = 0;
   
      while( ! ifDelim_ulCp(*ptrUL, delimUL) )
         ++ptrUL;
   
      tmpStr = (signed char *) ptrUL;
   
      while(*tmpStr != delimSC)
         ++tmpStr;
   
      return tmpStr - inStr;
   } /*lenStr_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun16: lenStrNull_ulCp
|   - finds the length of a string using ulong_ulCps
|     but also stops at null '\0'
| Input:
|   - inStr:
|     o c-string or string with deliminator to find length
|       of
|   - delimUL:
|     o deliminator (as ulong_ulCp (fun02)) at end of
|       string
|     o this will give an unused variable warning on
|       -DNOUL, but allows user to turn off ulCp
|   - delimSC:
|     o deliminator (as char) at end of string
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   lenStrNull_ulCp(
      signed char *inStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      signed char *tmpStr = inStr;

      delimSC = ~delimSC; 
      while(*tmpStr & delimSC)
         ++tmpStr;
   
      return tmpStr - inStr;
   } /*lenStrNull_ulCp*/
#else
   unsigned int
   lenStrNull_ulCp(
      signed char *inStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
      signed char *tmpStr = 0;
   
      while( ! ifNullDelim_ulCp(*ptrUL, delimUL) )
         ++ptrUL;
   
      tmpStr = (signed char *) ptrUL;

      delimSC = ~delimSC; 
      while(*tmpStr & delimSC)
         ++tmpStr;
      return tmpStr - inStr;
   } /*lenStrNull_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun17: endStr_ulCp
|   - finds the end of a c-string ('\0')
| Input:
|   - inStr:
|     o c-string or string to look for end in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   endStr_ulCp(
      signed char *inStr
   ){
      signed char *tmpStr = inStr;

      while(*tmpStr)
         ++tmpStr;
      return tmpStr - inStr;
   } /*endStr_ulCp*/
#else
   unsigned int
   endStr_ulCp(
      signed char *inStr
   ){
      ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
      signed char *tmpStr = 0;
   
      while( ! ifNull_ulCp(*ptrUL) )
         ++ptrUL;
   
      tmpStr = (signed char *) ptrUL;
   
      while( *tmpStr != '\0' )
         ++tmpStr;
   
      return tmpStr - inStr;
   } /*endStr_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun18: endLineUnix_ulCp
|   - finds the end of a c-string. This assumes that the
|     line ends in an '\0' or an '\n'
| Input:
|   - inStr:
|     o c-string or string to look for end in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   endLineUnix_ulCp(
      signed char *inStr
   ){
      signed char *tmpStr = inStr;
   
      while( *tmpStr & (~ '\n') )
         ++tmpStr;
      return tmpStr - inStr;
   } /*endLineUnix_ulCp*/
#else
   unsigned int
   endLineUnix_ulCp(
      signed char *inStr
   ){
      ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
      signed char *tmpStr = 0;
   
      while( ! ifLineUnix_ulCp(*ptrUL) )
         ++ptrUL;
   
      tmpStr = (signed char *) ptrUL;
   
      while( *tmpStr & (~ '\n') )
         ++tmpStr;
   
      return tmpStr - inStr;
   } /*endLineUnix_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun19: endLine_ulCp
|   - finds the end of a c-string (all OS's; looks for
|     '\0', '\n', and '\r')
|   - ingores all unused ascii characters (> 32, not '\t')
| Input:
|   - inStr:
|     o c-string or string to look for end in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   endLine_ulCp(
      signed char *inStr
   ){
      signed char *tmpStr = inStr;
   
      while(*tmpStr > '\r' || *tmpStr == '\t')
         ++tmpStr;
      return tmpStr - inStr;
   } /*endLine_ulCp*/
#else
   unsigned int
   endLine_ulCp(
      signed char *inStr
   ){
      ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
      signed char *tmpStr = 0;
   
      while( ! ifEndLine_ulCp(*ptrUL) )
         ++ptrUL;
   
      tmpStr = (signed char *) ptrUL;
   
      while(
            *tmpStr > '\r'  /*\r > \t > \n > \0*/
         || *tmpStr == '\t' /*so catch tab case*/
      ) ++tmpStr;
   
      return tmpStr - inStr;
   } /*endLine_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun20: endWhite_ulCp
|   - finds 1st white space ('\t', ' ', '\n', '\0', '\r')
| Input:
|   - inStr:
|     o c-string or string to look for whitespace in
| Output:
|   - Returns:
|     o number of characters in the string
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   endWhite_ulCp(
      signed char *inStr
   ){
      signed char *tmpStr = inStr;
   
      while( *tmpStr > 33)
         ++tmpStr;
      return tmpStr - inStr;
   } /*endWhite_ulCp*/
#else
   unsigned int
   endWhite_ulCp(
      signed char *inStr
   ){
      ulong_ulCp *ptrUL = (ulong_ulCp *) inStr;
      signed char *tmpStr = 0;
   
      while( ! ifWhite_ulCp(*ptrUL) )
         ++ptrUL;
   
      tmpStr = (signed char *) ptrUL;
   
      while( *tmpStr > 33)
         ++tmpStr;
   
      return tmpStr - inStr;
   } /*endWhite_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun21: eql_ulCp
|   - compares two strings until deliminator is found
| Input:
|   - qryStr:
|     o Pointer to query string
|   - refStr:
|     o Pointer to reference strin
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|     o this will give an unused variable warning on
|       -DNOUL, but allows user to turn off ulCp
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query > reference
|     o < 0 if query < reference
| Note:
|   - This will likely not be very good at comparing
|     short strings.
\-------------------------------------------------------*/
#ifdef NOUL
   signed long
   eql_ulCp(
      signed char *qryStr,
      signed char *refStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      while( *qryStr != delimSC )
      { /*Loop: find difference*/
         if(*qryStr != *refStr)
            break;
   
         ++qryStr;
         ++refStr;
      } /*Loop: find difference*/
   
      return *qryStr - *refStr;
   } /*eql_ulCp*/
#else
   signed long
   eql_ulCp(
      signed char *qryStr,
      signed char *refStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      ulong_ulCp *qryUL = (ulong_ulCp *) qryStr;
      ulong_ulCp *refUL = (ulong_ulCp *) refStr;
   
      while( ! ifDelim_ulCp(*qryUL, delimUL) )
      { /*Loop: Copy cpStr to dupStr*/
         if(*qryUL != *refUL)
            break;
   
         ++qryUL;
         ++refUL;
      } /*Loop: Copy cpStr to dupStr*/
   
      qryStr = (signed char *) qryUL;
      refStr = (signed char *) refUL;
   
      while( *qryStr != delimSC )
      { /*Loop: find difference*/
         if(*qryStr != *refStr)
            break;
   
         ++qryStr;
         ++refStr;
      } /*Loop: find difference*/
   
      return *qryStr - *refStr;
   } /*eql_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun22: eqlNull_ulCp
|   - compares two strings until null is found
| Input:
|   - qryStr:
|     o Pointer to query string
|   - refStr:
|     o Pointer to reference strin
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query > reference
|     o < 0 if query < reference
| Note:
|   - this will likely not be very good at comparing
|     short strings.
\-------------------------------------------------------*/
#ifdef NOUL
   signed long
   eqlNull_ulCp(
      signed char *qryStr,
      signed char *refStr
   ){
      while(*qryStr != (signed char) '\0' )
      { /*Loop: find difference*/
         if(*qryStr != *refStr)
            break;
   
         ++qryStr;
         ++refStr;
      } /*Loop: find difference*/
   
      return *qryStr - *refStr;
   } /*eqlNull_ulCp*/
#else
   signed long
   eqlNull_ulCp(
      signed char *qryStr,
      signed char *refStr
   ){
      ulong_ulCp *qryUL = (ulong_ulCp *) qryStr;
      ulong_ulCp *refUL = (ulong_ulCp *) refStr;
   
      while( ! ifNull_ulCp(*qryUL) )
      { /*Loop: Copy cpStr to dupStr*/
         if(*qryUL != *refUL)
            break;
   
         ++qryUL;
         ++refUL;
      } /*Loop: Copy cpStr to dupStr*/
   
      qryStr = (signed char *) qryUL;
      refStr = (signed char *) refUL;
   
      while( *qryStr != (signed char) '\0' )
      { /*Loop: find difference*/
         if(*qryStr != *refStr)
            break;
   
         ++qryStr;
         ++refStr;
      } /*Loop: find difference*/
   
      return *qryStr - *refStr;
   } /*eqlNull_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun23: eqlWhite_ulCp
|   - compares two strings until white space is found
| Input:
|   - qryStr:
|     o Pointer to query string
|   - refStr:
|     o Pointer to reference strin
| Output:
|   - Returns:
|     o 0 if strings are equal
|     o > 0 if query > reference
|     o < 0 if query < reference
| Note:
|   - this will likely not be very good at comparing
|     short strings.
\-------------------------------------------------------*/
#ifdef NOUL
   signed long
   eqlWhite_ulCp(
      signed char *qryStr,
      signed char *refStr
   ){
      while( *qryStr > 32 )
      { /*Loop: find difference*/
         if(*qryStr != *refStr)
            break;
   
         ++qryStr;
         ++refStr;
      } /*Loop: find difference*/
   
      if(*qryStr > 32)
         return *qryStr - *refStr;
      else if(*refStr > 32)
         return *qryStr - *refStr;
   
      return 0;
   } /*eqlWhite_ulCp*/
#else
   signed long
   eqlWhite_ulCp(
      signed char *qryStr,
      signed char *refStr
   ){
      ulong_ulCp *qryUL = (ulong_ulCp *) qryStr;
      ulong_ulCp *refUL = (ulong_ulCp *) refStr;
   
      while( ! ifWhite_ulCp(*qryUL) )
      { /*Loop: compare query to reference*/
         if(*qryUL != *refUL)
            break;
   
         ++qryUL;
         ++refUL;
      } /*Loop: compare query to reference*/
   
      qryStr = (signed char *) qryUL;
      refStr = (signed char *) refUL;
   
      while( *qryStr > 32 )
      { /*Loop: find difference*/
         if(*qryStr != *refStr)
            break;
   
         ++qryStr;
         ++refStr;
      } /*Loop: find difference*/
   
      if(*qryStr > 32)
         return *qryStr - *refStr;
      else if(*refStr > 32)
         return *qryStr - *refStr;
   
      return 0;
   } /*eqlWhite_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun24: rmWhite_ulCp
|   - removes white space from c-string
| Input:
|   - inStr:
|     o c-string to remove whites space
| Output:
|   - Modifies:
|     o  inStr to have white space removed
|   - Returns:
|     o number of characters in inStr
\-------------------------------------------------------*/
#ifdef NOUL
   unsigned int
   rmWhite_ulCp(
      signed char *inStr
   ){
      signed char *dupStr = inStr;
      signed char *cpStr = inStr;
   
      while(*cpStr != '\0')
      { /*Loop: remove white space*/
         if(*cpStr > 32)
            *dupStr++ = *cpStr++;
   
         else
         { /*Else: remove white space*/
            while(*cpStr < 33)
            { /*Loop: move past white space*/
               if(! *cpStr)
                  goto done_fun24;
               ++cpStr;
            } /*Loop: move past white space*/

            *dupStr++ = *cpStr++;
         } /*Else: remove white space*/
      } /*Loop: remove white space*/
   
      done_fun24:;
         *dupStr = '\0';
         return dupStr - inStr;
   } /*rmWhite_ulCp*/
#else
   unsigned int
   rmWhite_ulCp(
      signed char *inStr
   ){
      ulong_ulCp *cpUL = (ulong_ulCp *) inStr;
      ulong_ulCp *dupUL = (ulong_ulCp *) inStr;
   
      signed char *dupStr = inStr;
      signed char *cpStr = inStr;
   
      unsigned int uiChar = 0;
   
      while(*cpStr != '\0')
      { /*Loop: remove white space*/
         if(! ifWhite_ulCp(*cpUL) )
            *dupUL++ = *cpUL++;
   
         else
         { /*Else: have white space, manually copy*/
            cpStr = (signed char *) cpUL;
            dupStr = (signed char *) dupUL;
            
            for(
               uiChar = 0;
               uiChar < sizeof(ulong_ulCp);
               ++uiChar
            ){ /*Loop: remove white space*/
               if(*cpStr == '\0')
                  goto done_fun24;
               else if(*cpStr < 33)
                  ++cpStr;
               else
                  *dupStr++ = *cpStr++;
            } /*Loop: remove white space*/
   
            cpUL = (ulong_ulCp *) cpStr;
            dupUL = (ulong_ulCp *) dupStr;
         } /*Else: have white space, manually copy*/
      } /*Loop: remove white space*/
   
      done_fun24:;
         *dupStr = '\0';
         return dupStr - inStr;
   } /*rmWhite_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun25: swapDelim_ulCp
|   - swaps two strings until deliminator is found
| Input:
|   - firstStr:
|     o Pointer to string to first string to swap
|   - secStr:
|     o Pointer to second string to swap
|   - delimUL:
|     o delminator to end at (as long). Use makeULDelim
|       to build this deliminator
|     o this will give an unused variable warning on
|       -DNOUL, but allows user to turn off ulCp
|   - delimSC:
|     o delminator (as char) to stop copying at
| Output:
|   - Modifies:
|     o firstStr to have secStr string
|     o secStr to have firstStr string
| Note:
|   - This will likely not be very good at swapping short
|     strings.
\-------------------------------------------------------*/
#ifdef NOUL
   void
   swapDelim_ulCp(
      signed char *firstStr,
      signed char *secStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){
      while(
            *firstStr != delimSC
         && *secStr != delimSC
      ){ /*Loop: swap strings*/
         *firstStr ^= *secStr;
         *secStr ^= *firstStr;
         *firstStr++ ^= *secStr++;
      }  /*Loop: swap strings*/

      if(*firstStr != delimSC)
      { /*If: second string ended early*/
         *secStr++ = *firstStr;
         *firstStr++ = '\0';

         while(*firstStr != delimSC)
            *secStr++ = *firstStr++;
         *secStr = '\0';
      } /*If: second string ended early*/

      else if(*secStr != delimSC)
      { /*Else If: first string ended early*/
         *firstStr++ = *secStr;
         *secStr++ = '\0';

         while(*secStr != delimSC)
            *firstStr++ = *secStr++;
         *firstStr = '\0';
      } /*Else If: first string ended early*/

      else
      { /*Else: both strings were same length*/
         *firstStr = '\0';
         *secStr = '\0';
      } /*Else: both strings were same length*/
   } /*swapDelim_ulCp*/
#else
   void
   swapDelim_ulCp(
      signed char *firstStr,
      signed char *secStr,
      ulong_ulCp delimUL,
      signed char delimSC
   ){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
      ' Fun25: swapDelim_ulCp
      '   - swaps two strings until deliminator is found
      '   o fun25 sec01:
      '     - variable declarations
      '   o fun25 sec02:
      '     - swap until first deliminator
      '   o fun25 sec03:
      '     - if both strings end early
      '   o fun25 sec04:
      '     - if 1st string ends early, finsh swapping 2nd
      '   o fun25 sec05:
      '     - 2nd string ends early, finsh swapping 1st
      \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun25 Sec01:
      ^   - variable declarations
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
      ulong_ulCp *firstUL = (ulong_ulCp *) firstStr;
      ulong_ulCp *secUL = (ulong_ulCp *) secStr;
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun25 Sec02:
      ^   - swap until first deliminator
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
      while(
         ! (
               ifDelim_ulCp(*firstUL, delimUL)
             | ifDelim_ulCp(*secUL, delimUL) 
           )
      ){ /*Loop: Copy cpStr to dupStr*/
         *firstUL ^= *secUL;
         *secUL ^= *firstUL;
         *firstUL++ ^= *secUL++;
      } /*Loop: Copy cpStr to dupStr*/
   
      firstStr = (signed char *) firstUL;
      secStr = (signed char *) secUL;
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun25 Sec03:
      ^   - if both strings end early
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   

      if(
            ifDelim_ulCp(*firstUL, delimUL)
         && ifDelim_ulCp(*secUL, delimUL) 
      ){ /*If: both ended early*/
         while(
               *firstStr != delimSC
            || *secStr != delimSC
         ){ /*Loop: copy first string*/
            if(*secStr == delimSC)
            { /*If: second string ended first*/
               *secStr++ = *firstStr;
               *firstStr++ = '\0';
   
               while(*firstStr != delimSC)
                  *secStr++ = *firstStr++;
               *secStr = '\0';
               break;
            } /*If: second string ended first*/
   
            else if(*firstStr == delimSC)
            { /*Else If: first string ended first*/
               *firstStr++ = *secStr;
               *secStr++ = '\0';
   
               while(*secStr != delimSC)
                  *firstStr++ = *secStr++;
               *firstStr = '\0';
               break;
            } /*Else If: first string ended first*/
   
            else
            { /*Else: copying values*/
               *firstStr ^= *secStr;
               *secStr ^= *firstStr;
               *firstStr++ ^= *secStr++;
            } /*Else: copying values*/
         } /*Loop: copy first string*/
      }  /*If: both ended early*/
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun25 Sec04:
      ^   - if first string ends early, finsh swapping 2nd
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
      else if(ifDelim_ulCp(*firstUL, delimUL))
      { /*If: first string ended*/
   
         while(*firstStr != delimSC)
         { /*Loop: copy first string*/
            if(*secStr == delimSC)
               break;
   
            *firstStr ^= *secStr;
            *secStr ^= *firstStr;
            *firstStr++ ^= *secStr++;
         } /*Loop: copy first string*/
   
         if(*secStr != delimSC)
         { /*If: second string has more values*/
            *firstStr++ = *secStr;
            *secStr++ = '\0';
         } /*If: second string has more values*/
   
         firstUL = (ulong_ulCp *) firstStr;
         secUL = (ulong_ulCp *) secStr;
   
         while( ! ifDelim_ulCp(*secUL, delimUL) )
            *firstUL++ = *secUL++;
   
         firstStr = (signed char *) firstUL;
         secStr = (signed char *) secUL;
   
         while(*secStr != delimSC)
            *firstStr++ = *secStr++;
   
         *firstStr = '\0';
      } /*If: first string ended*/
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun25 Sec05:
      ^   - else 2nd string ends early, finsh swapping 1st
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
      else
      { /*Else: second string ended*/
   
         while(*secStr != delimSC)
         { /*Loop: copy first string*/
            if(*firstStr == delimSC)
               break;
   
            *firstStr ^= *secStr;
            *secStr ^= *firstStr;
            *firstStr++ ^= *secStr++;
         } /*Loop: copy first string*/
   
         if(*firstStr != delimSC)
         { /*If: first string has more values*/
            *secStr++ = *firstStr;
            *firstStr++ = '\0';
         } /*If: first string has more values*/
   
         firstUL = (ulong_ulCp *) firstStr;
         secUL = (ulong_ulCp *) secStr;
   
         while( ! ifDelim_ulCp(*firstUL, delimUL) )
            *secUL++ = *firstUL++;
   
         firstStr = (signed char *) firstUL;
         secStr = (signed char *) secUL;
   
         while(*firstStr != delimSC)
            *secStr++ = *firstStr++;
   
         *secStr = '\0';
      } /*Else: second string ended*/
   } /*swapDelim_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun26: swapNull_ulCp
|   - swaps two strings until null
| Input:
|   - firstStr:
|     o Pointer to string to first string to swap
|   - secStr:
|     o Pointer to second string to swap
| Output:
|   - Modifies:
|     o firstStr to have secStr string
|     o secStr to have firstStr string
| Note:
|   - This will likely not be very good at swapping short
|     strings.
\-------------------------------------------------------*/
#ifdef NOUL
   void
   swapNull_ulCp(
      signed char *firstStr,
      signed char *secStr
   ){
      while(*firstStr && *secStr)
      { /*Loop: swap strings*/
         *firstStr ^= *secStr;
         *secStr ^= *firstStr;
         *firstStr++ ^= *secStr++;
      }  /*Loop: swap strings*/

      if(*firstStr)
      { /*If: second string ended early*/
         *secStr++ = *firstStr;
         *firstStr++ = '\0';

         while(*firstStr)
            *secStr++ = *firstStr++;
         *secStr = '\0';
      } /*If: second string ended early*/

      else if(*secStr)
      { /*Else If: first string ended early*/
         *firstStr++ = *secStr;
         *secStr++ = '\0';

         while(*secStr)
            *firstStr++ = *secStr++;
         *firstStr = '\0';
      } /*Else If: first string ended early*/

      else
      { /*Else: both strings were same length*/
         *firstStr = '\0';
         *secStr = '\0';
      } /*Else: both strings were same length*/
   } /*swapNull_ulCp*/
#else
   void
   swapNull_ulCp(
      signed char *firstStr,
      signed char *secStr
   ){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
      ' Fun26: swapNull_ulCp
      '   - swaps two strings until deliminator is found
      '   o fun26 sec01:
      '     - variable declarations
      '   o fun26 sec02:
      '     - swap until first deliminator
      '   o fun26 sec03:
      '     - if both strings ended early
      '   o fun26 sec04:
      '     - if 1st string ends early, finsh swapping 2nd
      '   o fun26 sec05:
      '     - else 2nd string ends early, finsh 1st
      \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun26 Sec01:
      ^   - variable declarations
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
      ulong_ulCp *firstUL = (ulong_ulCp *) firstStr;
      ulong_ulCp *secUL = (ulong_ulCp *) secStr;
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun26 Sec02:
      ^   - swap until first deliminator
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
      while(
         ! (
                ifNull_ulCp(*firstUL)
              | ifNull_ulCp(*secUL)
           )
      ){ /*Loop: Copy cpStr to dupStr*/
         *firstUL ^= *secUL;
         *secUL ^= *firstUL;
         *firstUL++ ^= *secUL++;
      } /*Loop: Copy cpStr to dupStr*/
   
      firstStr = (signed char *) firstUL;
      secStr = (signed char *) secUL;
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun26 Sec03:
      ^   - if both strings ended early
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
      if(
           ifNull_ulCp(*firstUL)
        && ifNull_ulCp(*secUL)
      ){ /*If: both strings end early*/
         while(
               *firstStr != '\0'
            || *secStr != '\0'
         ){ /*Loop: copy first string*/
            if(*secStr == '\0')
            { /*If: second string ended first*/
               *secStr++ = *firstStr;
               *firstStr++ = '\0';
   
               while(*firstStr != '\0')
                  *secStr++ = *firstStr++;
               *secStr = '\0';
               break;
            } /*If: second string ended first*/
   
            else if(*firstStr == '\0')
            { /*Else If: first string ended first*/
               *firstStr++ = *secStr;
               *secStr++ = '\0';
   
               while(*secStr != '\0')
                  *firstStr++ = *secStr++;
               *firstStr = '\0';
               break;
            } /*Else If: first string ended first*/
   
            else
            { /*Else: copying values*/
               *firstStr ^= *secStr;
               *secStr ^= *firstStr;
               *firstStr++ ^= *secStr++;
            } /*Else: copying values*/
         } /*Loop: copy first string*/
      }  /*If: both strings end early*/
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun26 Sec04:
      ^   - if first string ends early, finsh swapping 2nd
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
        
      else if( ifNull_ulCp(*firstUL) )
      { /*Else If: first string ended*/
   
         while(*firstStr != '\0')
         { /*Loop: copy first string*/
            if(*secStr == '\0')
               break;
   
            *firstStr ^= *secStr;
            *secStr ^= *firstStr;
            *firstStr++ ^= *secStr++;
         } /*Loop: copy first string*/
   
         if(*secStr != '\0')
         { /*If: second string has more values*/
            *firstStr++ = *secStr;
            *secStr++ = '\0';
         } /*If: second string has more values*/
   
         firstUL = (ulong_ulCp *) firstStr;
         secUL = (ulong_ulCp *) secStr;
   
         while( ! ifNull_ulCp(*secUL) )
            *firstUL++ = *secUL++;
   
         firstStr = (signed char *) firstUL;
         secStr = (signed char *) secUL;
   
         while(*secStr != '\0')
            *firstStr++ = *secStr++;
   
         *firstStr = '\0';
      } /*Else If: first string ended*/
   
      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun26 Sec05:
      ^   - else 2nd string ends early, finsh swapping 1st
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
      else
      { /*Else: second string ended*/
   
         while(*secStr != '\0')
         { /*Loop: copy first string*/
            if(*firstStr == '\0')
               break;
   
            *firstStr ^= *secStr;
            *secStr ^= *firstStr;
            *firstStr++ ^= *secStr++;
         } /*Loop: copy first string*/
   
         if(*firstStr != '\0')
         { /*If: first string has more values*/
            *secStr++ = *firstStr;
            *firstStr++ = '\0';
         } /*If: first string has more values*/
   
         firstUL = (ulong_ulCp *) firstStr;
         secUL = (ulong_ulCp *) secStr;
   
         while( ! ifNull_ulCp(*firstUL) )
            *secUL++ = *firstUL++;
   
         firstStr = (signed char *) firstUL;
         secStr = (signed char *) secUL;
   
         while(*firstStr != '\0')
            *secStr++ = *firstStr++;
   
         *secStr = '\0';
      } /*Else: second string ended*/
   } /*swapNull_ulCp*/
#endif

/*-------------------------------------------------------\
| Fun27: shift_ulCp
|   - shifts a substring in a string up or down
| Input:
|   - inStr:
|     o c-string with substring to shift
|   - startSL:
|     o first character to shift (index 0)
|   - endSL:
|     o last character to shift (index 0)
|   - shiftSL:
|     o index to shift string at (index 0)
| Output:
|   - Modifies:
|     o inStr to have substring and startSL and endSL to
|       be at shiftSL
\-------------------------------------------------------*/
#ifdef NOUL
   void
   shift_ulCp(
      signed char *inStr, /*string with sub-string*/
      signed long startSL,/*first character to shift*/
      signed long endSL,  /*last character to shift*/
      signed long shiftSL /*were to start shift*/
   ){
      signed long lenSL = endSL - startSL + 1;
      signed long posSL = 0;

      if(startSL > shiftSL)
      { /*If: moveing sub-string back in string*/
         for( ; startSL <= endSL; ++startSL)
         { /*Loop: shift rotate item to target positon*/
            inStr[shiftSL] ^= inStr[startSL];
            inStr[startSL] ^= inStr[shiftSL];
            inStr[shiftSL] ^= inStr[startSL];
            ++shiftSL;

            for(
               posSL = shiftSL + lenSL;
               posSL < startSL;
               posSL += lenSL
            ){ /*Loop: shift old item correct position*/
               inStr[posSL] ^= inStr[startSL];
               inStr[startSL] ^= inStr[posSL];
               inStr[posSL] ^= inStr[startSL];
            }  /*Loop: shift old item correct position*/
         } /*Loop: shift rotate item to target positon*/
      } /*If: moveing sub-string back in string*/

      else
      { /*Else: moveing sub-string forward in string*/
         for( ; startSL <= endSL; ++startSL)
         { /*Loop: shift rotate item to target positon*/
            inStr[shiftSL] ^= inStr[startSL];
            inStr[startSL] ^= inStr[shiftSL];
            inStr[shiftSL] ^= inStr[startSL];
            ++shiftSL;

            for(
               posSL = startSL + lenSL;
               posSL < shiftSL;
               posSL += lenSL
            ){ /*Loop: shift old item correct position*/
               inStr[posSL] ^= inStr[startSL];
               inStr[startSL] ^= inStr[posSL];
               inStr[posSL] ^= inStr[startSL];
            }  /*Loop: shift old item correct position*/
         } /*Loop: shift rotate item to target positon*/
      } /*Else: moveing sub-string forward in string*/

      /*n * (distance / insertion_length)*/
   } /*shift_ulCp*/
#else
   void
   shift_ulCp(
      signed char *inStr, /*string with sub-string*/
      signed long startSL,/*first character to shift*/
      signed long endSL,  /*last character to shift*/
      signed long shiftSL /*were to start shift*/
   ){
      signed long lenSL = endSL - startSL + 1;
      signed long posSL = 0;

      unsigned long *inULPtr = 0;
      signed long lEndSL = endSL >> def_shiftULBy_ulCp;

      if(startSL > shiftSL)
      { /*If: moveing sub-string back in string*/
         inULPtr = (unsigned long *) &inStr[shiftSL];

         startSL >>= def_shiftULBy_ulCp;
         shiftSL >>= def_shiftULBy_ulCp;

         for( ; startSL <= lEndSL; ++startSL)
         { /*Loop: shift rotate item to target positon*/
            inULPtr[shiftSL] ^= inULPtr[startSL];
            inULPtr[startSL] ^= inULPtr[shiftSL];
            inULPtr[shiftSL] ^= inULPtr[startSL];
            ++shiftSL;

            for(
               posSL = shiftSL + lenSL;
               posSL < startSL;
               posSL += lenSL
            ){ /*Loop: shift old item correct position*/
               inULPtr[posSL] ^= inULPtr[startSL];
               inULPtr[startSL] ^= inULPtr[posSL];
               inULPtr[posSL] ^= inULPtr[startSL];
            }  /*Loop: shift old item correct position*/
         } /*Loop: shift rotate item to target positon*/

         startSL <<= def_shiftULBy_ulCp;
         shiftSL <<= def_shiftULBy_ulCp;

         for( ; startSL <= endSL; ++startSL)
         { /*Loop: shift rotate item to target positon*/
            inStr[shiftSL] ^= inStr[startSL];
            inStr[startSL] ^= inStr[shiftSL];
            inStr[shiftSL] ^= inStr[startSL];
            ++shiftSL;

            for(
               posSL = shiftSL + lenSL;
               posSL < startSL;
               posSL += lenSL
            ){ /*Loop: shift old item correct position*/
               inStr[posSL] ^= inStr[startSL];
               inStr[startSL] ^= inStr[posSL];
               inStr[posSL] ^= inStr[startSL];
            }  /*Loop: shift old item correct position*/
         } /*Loop: shift rotate item to target positon*/
      } /*If: moveing sub-string back in string*/

      else
      { /*Else: moveing sub-string forward in string*/
         inULPtr = (unsigned long *) &inStr[startSL];

         startSL >>= def_shiftULBy_ulCp;
         shiftSL >>= def_shiftULBy_ulCp;

         for( ; startSL <= lEndSL; ++startSL)
         { /*Loop: shift rotate item to target positon*/
            inULPtr[shiftSL] ^= inULPtr[startSL];
            inULPtr[startSL] ^= inULPtr[shiftSL];
            inULPtr[shiftSL] ^= inULPtr[startSL];
            ++shiftSL;

            for(
               posSL = startSL + lenSL;
               posSL < shiftSL;
               posSL += lenSL
            ){ /*Loop: shift old item correct position*/
               inULPtr[posSL] ^= inULPtr[startSL];
               inULPtr[startSL] ^= inULPtr[posSL];
               inULPtr[posSL] ^= inULPtr[startSL];
            }  /*Loop: shift old item correct position*/
         } /*Loop: shift rotate item to target positon*/

         startSL <<= def_shiftULBy_ulCp;
         shiftSL <<= def_shiftULBy_ulCp;

         for( ; startSL <= endSL; ++startSL)
         { /*Loop: shift rotate item to target positon*/
            inStr[shiftSL] ^= inStr[startSL];
            inStr[startSL] ^= inStr[shiftSL];
            inStr[shiftSL] ^= inStr[startSL];
            ++shiftSL;

            for(
               posSL = startSL + lenSL;
               posSL < shiftSL;
               posSL += lenSL
            ){ /*Loop: shift old item correct position*/
               inStr[posSL] ^= inStr[startSL];
               inStr[startSL] ^= inStr[posSL];
               inStr[posSL] ^= inStr[startSL];
            }  /*Loop: shift old item correct position*/
         } /*Loop: shift rotate item to target positon*/
      } /*Else: moveing sub-string forward in string*/

      /*n * (distance / insertion_length), the unsigned
      `   long trick will reduce time by 8x for long
      `   substrings
      */
   } /*shift_ulCp*/
#endif

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
