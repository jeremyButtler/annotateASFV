# Goal:

Give an idea on how to convert c-strings to base ten
  numbers with the base10str library.

Files: base10str.c base10str.h

# base10str

base10str is here to convert c-strings into nubmers.

Functions are in form of strTo<DataType>\_base10str. With
  <DataType> being the data type to convert to.

- Input:
  - c-string and a pointer to the data type to convert
  - `<data_type>` pointer to hold convert value
- Output:
  - returns number of characters converted to a number

| short hand | data Type      | function          |
|:----------:|:---------------|:------------------|
|     UL     | unsigned long  | strToUL\_base10str |
|     UI     | unsigned int   | strToUI\_base10str |
|     US     | unsigned short | strToUS\_base10str |
|     UC     | unsigned char  | strToUC\_base10str |
|     SL     | signed long    | strToSL\_base10str |
|     SI     | signed int     | strToSI\_base10str |
|     SS     | signed short   | strToSS\_base10str |
|     SC     | signed char    | strToSC\_base10str |
|     F      | float          | strToF\_base10str  |

Table: of base10str functions.

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "base10str.h"

int
main(
   int siArg,
   char *argAryStr[]
){
   signed long numberSL = 0;
   unsigned int indexUI = 0;

   if(siArg < 1)
   { /*If: nothing was input*/
      fprintf(stderr, "nothing input\n");
      return 1;
   } /*If: nothing was input*/

   indexUI =
      strToSL_base10str(
         (signed char *) argAryStr[1],
         &numberSL
      ); /*replace SL with any supported datatype you want
         `  to convert
         */

   if(argAryStr[1][indexUI])
   { /*If: failed to convert entire string*/
      fprintf(
         stderr,
         "string %s is non-numeric\n",
         argAryStr[1]
      );
      return 2;
   } /*If: failed to convert entire string*/

   else
      printf(
         "input: %s\noutput: %li\n",
         argAryStr[1],
         numberSL
     );

   return 0;
} /*main*/
```
