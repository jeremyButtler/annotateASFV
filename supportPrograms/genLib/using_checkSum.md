# Goal:

Learn how to use the functions in the checkSum library
  to calcualte a checksum.

Files: checkSum.h and checkSum.c

Dependencies: endin.h

# checkSum

Has functions to checksum checks. Currently I have the
  adler32 and crc32 checksum. For the crc32 functions,
  remember to set the initial crc32 value to -1.

For all checksums, you can add bytes to a checksum using
  the `byte` functions. However, at the end you must call
  the equivulent `finish` function to complete the
  checksum.

The main checkSum functions take an array of bytes
  (unsigned chars) and convert them to their checksum. For
  the non `us` functions the input is an array of bytes
  and the number of bytes to convert and the return value
  is a signed int
  (`checkSum = check_checkSum(bytesAry, length)`).

| Function            | Use                              |
|:--------------------|:---------------------------------|
| crc32\_checkSum     | find crc32 for a input string    |
| usAlder32\_checkSum | find adler32 for an input string |

Table: list of functions to get checksums from arrays of
  data.

The byte functions take in a byte and a checksum to add
  the byte to. This will be a signed integer
  (`checkByte_checkSum(byte, checkSum)`).

| Function                | Use                         |
|:------------------------|:----------------------------|
| crc32Byte\_checkSum     | add a byte to a crc32 sum   |
| adler32Byte\_checkSum   | add a byte to a adler32 sum |

Table: list of functions to add bytes to a current
  checkus (really on compiler to compile in fuction).
  When done adding bytes, you must call a finish function
  to finialize the checksum.

The Finish functions (call after added all bytes) take
  in the computed checksum (`idFinish_checkSum(checkSum)`
  and return the finialized checksum.

| Function                  | Use                       |
|:--------------------------|:--------------------------|
| crc32Finish\_checkSum     | do final step in a crc32  |
| adler32Finish\_checkSum   | do final step in adler32  |

Table: list of functions to finish off the checksum. These
  functions do not add bytes, but instead perform the
  final operation in computing the checksum.

Other functions include mkCrc32Tbl\_checkSum, which is not
  used, but can be used to auto-generate a crc32 table.
  Also, there is the crc32Mac\_checkSum macro in
  checkSum.c. The crc32Byte\_checkSum function is really a
  wrapper for crc32Mac\_checkSum.

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "checkSum.h"
#include "endin.h"

int
main(
){
   signed char *inputStr = (signed char *) "helow world!";
   signed int siChar = 0;

   /*for the crc32*/
   unsigned int crc32UI = -1; /*alawys needs to be -1*/

   /*for 32 bit addler calcuation*/
   signed int adler32SI = 0;
   signed int firstAdler32SI = 0;
   signed int secAdler32SI = 0;

   for(siChar = 0; inputStr[siChar] != '\0'; ++siChar)
   { /*Loop: find checksums*/
      crc32UI =
         crc32Byte_checkSum(inputStr[siChar], crc32UI); 
         /*getting crc32 calcuation*/

      adler32Byte_checkSum(
          inputStr[siChar],
          &firstAdler32SI,
          &secAdler32SI
      ); /*getting adler32 calcuation*/
   } /*Loop: find checksums*/

   /*finish off the crc32 calculation*/
   crc32UI = crc32Finish_checkSum(crc32UI);
   printf("crc32: %u\n", crc32UI);

   /*one round crc32 function*/
   printf(
      "crc32 auto %u\n",
      crc32_checkSum(inputStr,12)
   );


   /*finish off the 32 bit adler32 calcualtion*/
   adler32SI =
      adler32Finish_checkSum(firstAdler32SI,secAdler32SI);
   printf("adler32 32 bit: %i\n", adler32SI);

   /*one round adler32 function*/
   printf(
       "adler32 auto %i\n",
       adler32_checkSum(inputStr,12)
   );
} /*main*/
```
