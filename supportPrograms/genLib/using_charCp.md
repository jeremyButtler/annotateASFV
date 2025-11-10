# Goal:

Describe how to use the charCp library to work with
  c-strings.

Files: charCp.h and charCp.h

# charCp

The charCp library is not very good and the system
  libraries will generally be better. It is desinged
  to mainpulate c-strings byte by byte.

## length copy

cpLen\_charCp will copy x bytes from one string to
  another.

- Input:
  1. string to copy to
  2. string to copy from
  3. number of bytes to copy
- Output:
  1. Modifies input 1 to have copied bytes

```
signed char outStr[20];
signed char inStr[20] = (signed char *) "hellow world!";

cpLen_charCp(outStr, inStr, 13);
```

## copy till deliminator

You can copy until you reach a deimliminator using the
  `cpDelim_charCp()` function. Only be warned, this
  function will not stop for Nulls.

- Input:
  1. string to copy to
  2. string to copy from
  3. deliminator to stop at
- Output:
  1. Modfies input 1 to have copied string
  2. Returns number of bytes copied

```
signed char outStr[20];
signed char inStr[20] = (signed char *) "hellow world!";
unsigned int lenUI = 0;

lenUI = cpLen_charCp(outStr, inStr, ' ');
   /*copy till first space*/
```

## measure string lengths

You can measure the string lengths using
  the `lenStr_charCp()` and `lenStrNull_charCp()`
  functions. LenStrNull will measure untill a null or the
  input deliminator is found. lenStr will measure past
  nulls (stops only for input deliminator).

- Input:
  1. c-string to measure
  2. deliminator to stop at
- Output:
  1. Returns length of c-string till found deliminator

```
signed char inStr[20] = (signed char *) "hellow world!";
unsigned int lenUI = 0;

lenUI = lenStr_charCp(inStr, ' ');
   /*returns 7 (length of hellow)*/

lenUI = lenStr_charCp(inStr, 'x');
   /*goes out of bounds, so no idea on return; */

lenUI = lenStrNull_charCp(inStr, 'x');
   /*returns 13 (length of "hellow world!")*/
```

## compare strings

You can compare strings using the `eql_charCp()` and
  `match_charCp()` functions. Eql will compare two strings
  untill a deliminator is found or a difference is found.
  Match, is like eql, but it will assume both strings are
  equal if only the query ends early.

- Input eql:
   1. first c-string to compare
   2. second c-string to compare
   3. deliminator to stop at
- Output:
   1. 0 for equal
   2. < 0 if input 1 is less than input 2
   3. > 0 if input 1 is greater than input 2

- Input match:
   1. query c-string to compare
   2. reference c-string to compare
   3. deliminator to stop at
   4. signed int pionter to hold index of when comparision
      stopped
- Output:
   1. 0 for equal
   2. < 0 if input 1 is less than input 2
   3. > 0 if input 1 is greater than input 2


```
signed char oneStr[20] = (signed char *) "hellow world!";
signed char twoStr[20] =
   (signed char *) "hellow world! by";
signed int lastCharSI = 0;

if(! eql_charCp(oneStr, twoStr, ' ') )
   printf("equal for space\n");
   /*in this case will be true since both strings end at
   `  the same point
   */

if(! eql_charCp(oneStr, twoStr, '\0') )
   printf("equal for space\n");
   /*this will be false; since twoStr has "by"*/

if(! match_charCp(oneStr, twoStr, &lastCharSI, '\0') )
   printf("match for null\n");
   /*this will be true; since the query is shorter then
   `  the referece
   `lastCharSI will be set to 13 (end of query)
   */
```
