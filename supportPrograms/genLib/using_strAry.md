# Goal:

Describe how to use strAry to store and get multiple
  c-strings from a c-string string array.

Files: strAry.h and strAry.c

Dependencies: ulCp

# strAry

The strAry library stores multiple c-string's in a signed
  character array. Each c-string can use a maximum of 63
  bytes. In all cases c-strings take up the 64 bytes
  (63 bytes of data and 1 byte for null).

There is no structure, since the modified/returned values
  are a signed integer and a signed character array.

## making and freeing c-string arrays

You can make a c-string array using the `mk_strAry()`
  function to allocate the initial memory.

- Input:
  - number of c-strings to store. you can reallocate later
- Output:
  - Returns:
    - heap allocated c-string array to use
    - 0 for memory errors

You can then free the c-string array using the `free()`
  function in the standard library (like with malloc and
  calloc).

```
signed char *strAry = mk_strAry(10);
   /*stores 10 c-strings*/
if(! strAry)
  /*had memory error*/

if(strAry)
   free(strAry);
strAry = 0;
```

## adding strings, getting strings, and increasing memory

The c-string array made by `mk_strAry()` is a single
  dimensional array that is split into 64 byte chunks.
  This means adding c-strings, getting c-strings, or
  increasing the memory is often complicated. So, use
  the base functions.

You can add a c-string to a c-string array using the
  `add_strAry()` function.

- Input:
  1. c-string to add in
  2. c-string array (signed char \*) to add input 1 to
  3. index to add input 1 to (make sure you have room)

You can get a pointer to a string using the `get_strAry()`
  function. This is here because the array is split into
  64 byte chunks. So, index 1 is actually 65.

- Input:
  1. c-string array (signed char \*) to get pointer from
  2. index of c-string to get
- Output:
  - Returns: pointer to c-string at index of input 2

You can increase the memory of a c-string array by using
  the `realloc_strAry()` function. This does not call
  malloc, so, make sure to use `mk_strAry()` first. This
  is really just a wrapper around the realloc() function.

- Input:
  1. c-string array to increase memory of
  2. number of c-strings for the arrays to add
- Output:
  - Returns:
    - pointer to resized c-string array
      - You must free the old array
    - 0 for memory errors

**Example**

```
signed char *swapStr = 0;
signed char *strAry = mk_strAry(2);

add_strAry(strAry, (signed char *) "hellow world!", 0);
add_strAry(strAry, (signed char *) "by world!", 1);

swapStr = realloc_strAry(strAry, 3);
if(! swapStr)
   /*deal with memory errors*/
else
   strAry = swapStr;

add_strAry(strAry, (signed char *) "by world!", 1);

printf("string at index 1 is %s\n", get_strAry(strAry,1));
   /*this would get the pointer for the string at index 1
   `  "by world!"
   */
```

## sorting

You can sort your c-string array with the sort function
  (uses shell sort) or by adding strings in a sorted
  order (using a binrary search, but still need to move
  strings around, so O(n^2)).

### sorting during addition

To add strings in a sorted order you must use
  the `addSort_strAry()` function. You can not use
  `add_strAry()`, unless you sort the strings before
  calling `addSort_strAry()`. This will also mess up
  your index.

- WARNINGs:
  - this does not rellocate memory, so make sure you have
    room for one more string
  - is this is a O(n^2) time do to shifting strings around

- Input:
  1. c-string to add in
  2. c-string array to add input 1 to
  3. number elements in the c-string array (input 2)
- Output:
  - Returns: index of added string

**Example**

```
signed char *strAry = mk_strAry(2);

addSort_strAry(strAry,(signed char *) "hellow world!",0);
addSort_strAry(strAry, (signed char *) "by world!", 1);

if(strAry)
   free(strAry);
strAry = 0;
```

### sorting when finished

The other option is to add strings and then sort them
  later with the `sort_strAry()` function.

- Input:
  1. c-string to sort
  2. length (index 1) of c-string to sort

You keep an unsigned int array in sync with the strings
  using the `sortSync_strAry()` function.

- Input `sortSync_str_ptrAry()`:
  1. c-string to sort
  2. unsigned int array to keep in sync
  3. length (index 1) of c-string to sort

You keep an unsigned short array in sync with the strings
  using the `sortSyncUS_strAry()` function.

- Input `sortSync_str_ptrAry()`:
  1. c-string to sort
  2. unsigned short array to keep in sync
  3. length (index 1) of c-string to sort


**Example; regular sorting**

```
signed char *swapStr = 0;
signed char *strAry = mk_strAry(2);

add_strAry(strAry, (signed char *) "hellow world!", 0);
add_strAry(strAry, (signed char *) "by world!", 1);

sort_strAry(strAry, 2);

if(strAry)
   free(strAry);
strAry = 0;
```

**Example; keeping an array in sync sorting**

```
signed char *swapStr = 0;
signed char *strAry = mk_strAry(2);
unsigned int *uiAry = {0, 1};

add_strAry(strAry, (signed char *) "hellow world!", 0);
add_strAry(strAry, (signed char *) "by world!", 1);

sortSync_strAry(strAry, uiAry, 2);

if(strAry)
   free(strAry);
strAry = 0;
```

## finding c-strings

You can find c-strings using the `find_strAry()` and
  the `findNoSort_strAry()` functions. Both functions take
  the same input and return the same output.

The difference is that `find_strAry()` does a binary
  search and needs a sorted array,
  while `findNoSort_strAry()` starts from the start
  and moves till the string is found or the end of the
  array.

- Input (both):
  1. c-string array to search
  2. c-string to find in the c-string array
  3. length (index 1) of c-string to sort
- Output:
  - Returns:
    - index of found c-string
    - -1 if could not find the c-string

**Example**

```
signed char *swapStr = 0;
signed char *strAry = mk_strAry(2);

add_strAry(strAry, (signed char *) "hellow world!", 0);
add_strAry(strAry, (signed char *) "by world!", 1);

sort_strAry(strAry, 2);

printf(
   "by world! is at index %li",
   find_strAry(strAry, (signed char *) "by world!", 2)
); /*returns 1 (index)*/

printf(
   "this returns index -1 (string not in array) %li",
   find_strAry(strAry, (signed char *) "not orld!", 2)
); /*returns 1 (index)*/

if(strAry)
   free(strAry);
strAry = 0;
```

## misc functions

### string comparison

You can compare two strings in a str\_ptrAry structure
  using the `cmp_strAry()` to compare a single string to
  a string in a string array. You can use
  the `cmpIndex_strAry()` function to compare two strings
  in a string array.

- Input (`cmp_strAry()`):
  1. c-string to compare (query)
  2. c-string array with string to compare to
  3. index of string to compare to (reference)
- Output:
  - Returns:
    - 0 if equal
    - > 0 if query (input 1) is greater
    - < 0 if reference (input 3) is greater

- Input (`cmpIndex_strAry()`):
  1. c-string array with string to compare to
  1. index of first string to compare (queyr)
  3. index of string to compare to (reference)
- Output:
  - Returns:
    - 0 if equal
    - > 0 if query (input 2) is greater
    - < 0 if reference (input 3) is greater

```
signed char *swapStr = 0;
signed char *strAry = mk_strAry(2);

add_strAry(strAry, (signed char *) "hellow world!", 0);
add_strAry(strAry, (signed char *) "hellow world!", 1);

/*comparing two index's*/
if(! cmpIndex_strAry(strAry, 0, 1) )
   add_strAry(strAry, (signed char *) "by world!", 1);
   /*means index 0 and 1 are the same*/

if(! cmp_strAry((signed char *) "by world!", strAry, 1) )
   printf("index 1 is: by world!\n");

if(strAry)
   free(strAry);
strAry = 0;
```

### string swapping

You can swap two strings in a c-string array using
  the `swap_strAry()` function.

- Input:
  1. c-string array to swap elements in
  2. index of first string to swap
  3. index of second string to swap

```
signed char *swapStr = 0;
signed char *strAry = mk_strAry(2);

add_strAry(strAry, (signed char *) "hellow world!", 0);
add_strAry(strAry, (signed char *) "by world!", 1);

swap_strAry(strAry, 0, 1);
   /*move "by world!" to "hellow world!"*/

if(strAry)
   free(strAry);
strAry = 0;
```
