# Goal:

Describe how to convert numbers to c-strings with the
  numToStr library.

Files: numToStr.h and numToStr.h

# numToStr:

Converts unsigned numbers (numToStr) or doubles
  (double\_numToStr) to c-strings. For doubles, negatives
  values are allowed.

For signed integers use the `signed_numToStr()` function.

If you want an unsigned integer to be copied bacwards,
  use the `backwards_numToStr()` function.

- Input (all except double\_numToStr):
  1. c-string to add number to
  2. number to convert
- Output (all functions)
  1. Modifies input 1 to have converted number
  2. returns length of input 1

- Input double\_numToStr:
  1. c-string to add number to
  2. number to convert
  3. maximum number of decimal digits
- Output (all functions)
  1. Modifies input 1 to have converted number, rounded
     up from input 3
  2. returns length of input 1

**Example**

```
signed char outStr[20];
unsigned int digitsConvertedUI = 0;

digitsConvertedUI = numToStr(outStr, 100);
   /*digitsConvertedUI will be 3, and outStr will
   `  be "100\0"
   */
digitsConvertedUI = numToStr(outStr, -1);
   /*digitsConvertedUI will be 9 (max digits for integer),
   `  and outStr will have the maximum value for an
   `  integer
   */


digitsConvertedUI = signed_numToStr(outStr, -100);
   /*digitsConvertedUI will be 4, and outStr will
   `  be "-100\0"
   */


digitsConvertedUI = double_numToStr(outStr, 0.15, 2);
   /*digitsConvertedUI will be 4, and outStr will
   `  be "0.15\0"
   */
digitsConvertedUI = double_numToStr(outStr, 0.15, 1);
   /*digitsConvertedUI will be 4, and outStr will
   `  be "0.2\0"
   */
```
