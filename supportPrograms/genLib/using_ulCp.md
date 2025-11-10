# Goal:

Tell how to copy, get string lengths, compare, and
  mainpulate strings with the ulCp library.

Files: ulCp.h, ulCp.c, and 64bit.h

# ulCp: overview

ulCp is short for unsigned long copy. It is like charCp,
  except that it uses unsigned long. As a general rule it
  is slower thang strcpy, which uses SIMD. However, it
  can be faster than charCp or other byte by byte
  methods (at least on 64 bit) when strings are longer
  then 16 bytes. It also allows for you to use your own
  deliminators for some functions.

Use mkDelim\_ulCp to make a deliminator and store it in
  an ulong\_ulCp dataType.

There are a lot of functions in ulCp, so I am not going
  to go over them all. They fall into deliminator making
  functions, detection functions, copying functions,
  string length functions, string comparision functions,
  white space clean up, and swap functions.

Note: in swap functions, the `shift_ulCp()` was never
  tested, so likely does not work.

# ulCp Macros

These macros are designed to check a singe unsigned long
  length of bytes in a string for a particular search
  pattern.

The one exception here is mkDelim, which makes the
  deliminator.

## deliminator making function

The deliminator marks the end of the string for functions
  that take a deliminator.

You can make a deliminator with the `mkDelim_ulCp()`
  function.

- Input:
  1. character to use as the deliminator
- Output:
  - Returns:
    - ulong\_ulCp datatype with the deliminator to use

```
ulong_ulCp delimUL = mkDelim_ulCp('A');
```

## checking for deliminators

These functions (macros) check to see if deliminators are
  in the next 8 bytes (less if unsigned long is less then
  8 bytes) in a string.

### deliminator check

The `ifDelim_ulCp()` looks to see if the deliminator
  from `mkDelim_ulCp()` is the next few bytes of the
  string. The string should be in the from of an
  ulong\_ulCp datatype.

- Input:
  - c-string as a ulong\_ulCp pointer to check
  - deliminator (as ulong\_ulCp data type) to look for
- Output:
  - 0 if deliminator not found
  - not 0 if deliminator was found

```
signed char *checkStr = "hellow world!";
ulong_ulCp *ulStr = (ulong_ulCp *) checkStr;
ulong_ulCp delimUL = mkDelim_ulCp(' ');

if( ifDelim_ulCp(ulStr, delimUL) )
   printf("space is in first 8 chars of %s", checkStr);
   /*this will fire*/
```

### deliminator and null check

The `ifNullDelim_ulCp()` looks to see if null or the
  deliminator from `mkDelim_ulCp()` is the next few bytes
  of the string. The string should be in the from of an
  ulong\_ulCp datatype.

- Input:
  - c-string as a ulong\_ulCp pointer to check
  - deliminator (as ulong\_ulCp data type) to look for
- Output:
  - 0 if deliminator not found
  - not 0 if deliminator was found

```
signed char *checkStr = "hellow world!";
ulong_ulCp *ulStr = (ulong_ulCp *) checkStr;
ulong_ulCp delimUL = mkDelim_ulCp(' ');

if( ifNullDelim_ulCp(ulStr, delimUL) )
   printf("space or null in first part of %s", checkStr);
   /*this will fire*/
```

### null check

The `ifNull_ulCp()` looks to see if null is in the next
  few bytes of the string. The string should be in the
  from of an ulong\_ulCp datatype.

- Input:
  - c-string as a ulong\_ulCp pointer to check
- Output:
  - 0 if deliminator not found
  - not 0 if deliminator was found

```
signed char checkStr[] = "hellow world!";
ulong_ulCp *ulStr = (ulong_ulCp *) checkStr;

if( ifNull_ulCp(ulStr) )
   printf("null in first part of %s", checkStr);
   /*this will not fire*/

else
{ /*Else: make case that will fire*/
   checkStr[7] = 0;
   
   if( ifNull_ulCp(ulStr) )
      printf("null in first part of %s", checkStr);
      /*this will fire*/
} /*Else: make case that will fire*/
```

### check end of lines

You can check for line endings using
  the `ifLineUnix_ulCp()` or the `ifEndLine_ulCp()`
  fuctions. The string should be in the form of an
  ulong\_ulCp datatype.

The `ifLineUnix_ulCp()` function looks for `\n` or
  nulls (`\0`).

The `ifEndLine_ulCp()` function looks for `\n`, `\r`, or
  nulls (`\0`). So, should cover all OS's.

- Input (both):
  - c-string as a ulong\_ulCp pointer to check
- Output:
  - 0 if deliminator not found
  - not 0 if deliminator was found

```
signed char checkStr[] = "hello\n";
ulong_ulCp *ulStr = (ulong_ulCp *) checkStr;

if( ifEndLine_ulCp(ulStr) )
   printf("%s line ends in next 8 bytes", checkStr);
   /*this will fire*/
```

### check for white space

If you want to look for white space (ascii characters <
  33) you can use `ifWhite_ulCp()`. The string should be
  in the form of an ulong\_ulCp data type.

- Input (both):
  - c-string as a ulong\_ulCp pointer to check
- Output:
  - 0 if white space was not found
  - not 0 if white space was found

```
signed char checkStr[] = "hello world!";
ulong_ulCp *ulStr = (ulong_ulCp *) checkStr;

if( ifWhite_ulCp(ulStr) )
   printf("%s has white space in next 8 bytes", checkStr);
   /*this will fire*/
```

# ulCp functions

These functions are designed to copy, get length, edit,
  or compare strings. The call the ulCp macros for
  the comparisions.

## string copying

These functions are used to copy strings from one string
  to another.

### copy string by length

The `cpLen_ulCp()` function copies a string by input
  length, using unsigned long.

- Input:
  1. c-string to copy to
  2. c-string to copy from
  3. number of bytes to copy
- Output:
  1. fills input 1 with characters from input 2


```
signed char outStr[32];
signed char *copyStr = "hello world!";

cpLen_ulCp(outStr, copyStr, 10);
printf("%s", outStr);
   /*this will print the first 10 characters, which
   `  are "hello worl"
   */
```

### copy string by deliminator

The `cpDelim_ulCp()` function copies a string until the
  input deliminator is found.

- Input:
  1. c-string to copy to
  2. c-string to copy from
  3. ulong\_ulCp deliminator
  4. signed char deliminator
- Output:
  - Modifies: input 1 to have string in input 2
  - Returns: number of characters copied
- Warning: will ignore nulls


```
signed char outStr[32];
signed char *copyStr = "hello world!";
ulong_ulCp delimUL = mkDelim_ulCp('o');
unsigned int lenUI = 0;

lenUI = cpDelim_ulCp(outStr, copyStr, delimUL, 'w');
   /*lenUI will be set to 7 (bytes copied)*/
printf("%s", outStr);
   /*this will print "hello w" (stops at o)*/
```

### copy entire string

The `cpStr_ulCp()` function copies a string.

- Input:
  1. c-string to copy to
  2. c-string to copy from
- Output:
  - Modifies: input 1 to have string in input 2
  - Returns: number of characters copied

```
signed char outStr[32];
signed char *copyStr = "hello world!";
unsigned int lenUI = 0;

lenUI = cpStr_ulCp(outStr, copyStr);
   /*lenUI will be set to 7 (bytes copied)*/
printf("%s", outStr);
   /*this will print "hello world!"*/
```

### copy a line until line break

The `cpLineUnix_ulCp()` function copies a string untill
  a Unix line break (`\n`) or null (`\0`) is found.

You can use the `cpLine_ulCp()` function to copy a line
  untill a line break (`\n` or `\r`) or null (`\0`) is
  found.

- Input (cpLineUnix and cpLine):
  1. c-string to copy to
  2. c-string to copy from
- Output:
  - Modifies: input 1 to have string in input 2
  - Returns: number of characters copied

```
signed char outStr[32];
signed char *copyStr = "hello world!\n";
unsigned int lenUI = 0;

lenUI = cpLineUnix_ulCp(outStr, copyStr);
printf("%s", outStr);
   /*this will print "hello world!"*/
```

### copy a line till white space found

The `cpWhite_ulCp()` function copies a string untill
  a white space
  (ascii < 33 [space, tab, line break, null]) is found.

- Input (cpLineUnix and cpLine):
  1. c-string to copy to
  2. c-string to copy from
- Output:
  - Modifies: input 1 to have string in input 2
  - Returns: number of characters copied

```
signed char outStr[32];
signed char *copyStr = "hello world!\n";
unsigned int lenUI = 0;

lenUI = cpWhite_ulCp(outStr, copyStr);
printf("%s", outStr);
   /*this will print "hello"*/
```

## finding string length

You can find the string lengths using functions in ulCp.

### finding length by deliminators

You can use the `lenStr_ulCp()` function to find the
  string length for a deliminator. The warning is that
  this function does not check for nulls.

To check for nulls and deliminators, use
  the `lenStrNull_ulCp()` function.

- Input (lenStr and lenStrNull):
  1. c-string to get length of
  2. ulong\_ulCp datatype with the deliminator
  3. signed char with the ascii deliminator
- Output:
  1. Returns length of c-string till the deliminator

```
signed char *inStr = "hello world!";
ulong_ulCp delimUL = mkDelim_ulCp('r');
unsigned int lenUI = 0;

lenUI = lenStr_ulCp(outStr, copyStr, delimUL, 'r');
printf("string is %u long\n", lenUI);
   /*this will print 9 (measures till first r)"*/

/*this varation will stop at nulls if they are present*/
lenUI = lenStrNull_ulCp(outStr, copyStr, delimUL, 'r');
printf("string is %u long\n", lenUI);
   /*this will print 9 (measures till first r)"*/
```

### find c-string length

You can find the c-string length using the `endStr_ulCp()`
  function. This stops at a null.

- Input:
  1. c-string to measure
- Output:
  1. Returns length of c-string

```
signed char *inStr = "hello world!";
unsigned int lenUI = 0;

lenUI = endStr_ulCp(outStr, copyStr);
printf("string is %u long\n", lenUI);
   /*this will print 12 (length of "hello world!")*/
```

### find length of one line in c-string

You can find the length of one line (till line break)
  using the `endLineUnix_ulCp()` or `endLine_UlCp()`
  functions. The endLineUnix function looks for a new
  line (`\n`) or null (`\0`), while endLine looks for
  the line break (`\n`), carriage return (`\r`), or a
  null (`\0`). So, endLineUnix is for Unix, while
  endLine should work for most OS's.

- Input:
  1. c-string to measure
- Output:
  1. Returns length of c-string

```
signed char *inStr = "hello world!\r\nextra stuff";
unsigned int lenUI = 0;

lenUI = endLineUnix_ulCp(outStr, copyStr);
printf("string is %u long\n", lenUI);
   /*this will print 13 (length of "hello world!\r")*/

lenUI = endLine_ulCp(outStr, copyStr);
printf("string is %u long\n", lenUI);
   /*this will print 12 (length of "hello world!")*/
```

### find first white space in c-string

You can find the first white space in a string using
  the `endWhite_ulCp()` function. This looks for spaces,
  tabs, newlines, and nulls. It will also search for
  any ascii character less then 33 (non-used).

- Input:
  1. c-string to measure
- Output:
  1. Returns length of c-string

```
signed char *inStr = "hello world!";
unsigned int lenUI = 0;

lenUI = endWhite_ulCp(outStr, copyStr);
printf("string is %u long\n", lenUI);
   /*this will print 5 (length of "hello")*/
```

## comparsion functions

The comparison functions are here to compare two strings
  to each other. The function type determines what
  deliminator is used.

### variable deliminator comparisions

You can use the `eql_ulCp()` function to compare two
  strings untill the input deliminator is found. However,
  this will not stop at nulls.

- Input:
  1. query c-string to compare
  2. reference c-string to compare
  3. ulong\_ulCp with deliminator to stop at
  4. signed char with deliminator to stop at
- Output:
  - Returns: 0 if strings are equal
  - Returns: < 0 if query is less then reference
  - Returns: > 0 if query is greater then reference

```
signed char *queryStr = "hello world!";
signed char *referenceStr = "hello yes!";
ulong_ulCp delimUL = mkDelim_ulCp(' ');
signed long resultSL = 0;

resultSL = eql_ulCp(queryStr, referenceStr, delimUL, ' ');
   /*this will return 0, since both strings are matches
   `  until the first space
   */
```

### comparing until null

You can compare c-strings untill a null is found with
  the `eqlNull_ulCp()` function.

- Input:
  1. query c-string to compare
  2. reference c-string to compare
- Output:
  - Returns: 0 if strings are equal
  - Returns: < 0 if query is less then reference
  - Returns: > 0 if query is greater then reference

```
signed char *queryStr = "hello world!";
signed char *referenceStr = "hello world!";
signed long resultSL = 0;

resultSL = eqlNull_ulCp(queryStr, referenceStr);
   /*this will return 0, since both strings are the same*/
```

### comparing until white space

You can use the `eqlWhite_ulCp()` function to compare two
  strings untill any white space (ascii less then 33) is
  found [ex: space, tab, carriage return, new line, null].

- Input:
  1. query c-string to compare
  2. reference c-string to compare
- Output:
  - Returns: 0 if strings are equal
  - Returns: < 0 if query is less then reference
  - Returns: > 0 if query is greater then reference

```
signed char *queryStr = "hello world!";
signed char *referenceStr = "hello yes!";
signed long resultSL = 0;

resultSL = eqlWhite_ulCp(queryStr, referenceStr);
   /*this will return 0, since both strings are matches
   `  until the first space
   */
```

## string clean up

Currently the only string clean up function is
  the `rmWhite_ulCp()` function. This removes white space
  from the input string. It uses unsigned longs to move
  past regions without whites space and byte by byte
  editing to remove white space.

- Input:
  1. c-string to remove white space from
- Output:
  1. Modifies: input 1 to have no white space
  2. Returns: length of cleaned up string

```
signed char inStr[32] = "hello world!\r\n\tend";
unsigned int lenUI = 0;

lenUI = rmWhite_ulCp(inStr);
printf("length: %u; string: %s\n", lenUI, inStr);
   /*output should be
   `   "length: 15; string: helloworld!end\n"
   */
```

## swaping functions

There are a couple functions to swap c-strings.

### swaping by deliminator

You can swap c-strings using deliminators using
  the `swapDelim_ulCp()` function. It will not stop at
  null, unless the deliminator is null.

- Input:
  1. first c-string to swap
  2. second c-string to swap
  3. ulong\_ulCp with delimator to stop at
  4. delimSC signed char with delimator to stop at
- Output:
  1. Modifies: first c-string (input 1) to have the
     contents of second c-string (input 2)
  2. Modifies: second c-string (input 2) to have the
     contents of first c-string (input 1)

```
signed char firstStr[32] = "hello world!";
signed char secondStr[32] = "by world!";
ulong_ulCp delimUL = mkDelim_ulCp(' ');

swapDelim_ulCp(firstStr, secondStr, delimUL, ' ');
   /*swap strings until the space is found, so
   `   firstStr is "by"
   `   secondStr is "hello"
   */
```

### swaping full strings

You can swap entire strings using the `swapNull_cpDelim()`
  function.

- Input:
  1. first c-string to swap
  2. second c-string to swap
- Output:
  1. Modifies: first c-string (input 1) to have the
     contents of second c-string (input 2)
  2. Modifies: second c-string (input 2) to have the
     contents of first c-string (input 1)

```
signed char firstStr[32] = "hello world!";
signed char secondStr[32] = "by world!";

swapNull_ulCp(firstStr, secondStr);
   /*swap full strings
   `   firstStr is "by world!"
   `   secondStr is "hello world!"
   */
```

# final notes

## ulCp: warnings

I have no idea if this will work on big endin systems. Not
  a big deal since most user OSs are little endin. I have
  tested one or two functions using qemu to emulate a
  mk68040 processor on a debian image I found online
  ([https://people.debian.org/~gio/dqib/](
   https://people.debian.org/~gio/dqib/)
  ). It seems to work, but then I also had mixed results
  with a browser emulator.

Some CPUs or OS's (not sure which) will require integers
  to be 32 bit aligned, which means casting characters
  to longs will not work. If you are on an OS were
  compatibility is an issue, then use `-DNOUL` to do
  byte by byte copying. You will get unused variable
  warnings for functions that use a long deliminator and
  a character deliminator.

## ulCp: some details

I had to define my own data type for ulCp because plan9
  treats longs as 32 bits. So, their is a special flag
  (`-DPLAN9_64`) which will use `long long`
  as `ulong_ulCp` to get 64 bit dataTypes (only use on
  64 bit systems).

I have defined several deliminators at the top of the
  ulCp.h file. To see a list of functions vew the table
  of contents (Start Of File/SOF) at the start of ulCp.c
  or ulCp.h. I define 64 bit values and then use casting
  to resize it to 32 bit or whatever size is needed.

| deliminator     | defined macro          |
|:----------------|:-----------------------|
| null            | def\_null\_ulCp (or 0) |
| new line        | def\_newline\_ulCp     |
| carriage return | def\_carriage\_ulCp    |
| tab             | def\_tab\_ulCp         |
| space           | def\_space\_ulCp       |
| comma           | def\_comma\_ulCp       |

Table: of deliminators in ulCp.h

There are also several long check macros to see if a value
  is in a long (fun03 to fun08 in ulCp.h). Using these
  with `-DNOUL` will be very inefficent. They really on
  a four to eight bytes being stored in longs.

| deliminator     | function           |
|:----------------|:------------------:|
| user defined    | ifDelim\_ulCp      |
| null + user     | ifDelimNull\_ulCp  |
| null ('\\0')    | ifNull\_ulCp       |
| newline + null  | ifLineUnix\_ulCp   |
| \\n + \\r + \\0 | ifEndLine\_ulCp    |
| white space     | ifWhite\_ulCp      |

Table: of macros that check a long for certain charaters.
  for user or user defined deliminators, use a pre-defined
  deliminator or `ulong_ulCp delim = mkDelim_ulCp(value)`
  to make a deliminator. By white space I mean any
  asccii value less then 33 (this includes space, tab,
  carriage return, newline, and null).

- Flags:
  - `-DNOUL`: do byte by byte coping, slower, but always
    works byte copying
  - `-DPLAN9_64` for plan9 longs are 32bit, so you need ot
    use `long long` on 64 bit plan9 OS's

Finally, this has a slew of functions, check the header
  of the .c or .h file for a list.
