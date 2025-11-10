# Goal:

Desribe how to use ptrAry to handle dynamic c-string
  arrays (a c-string pionter array).

Files: ptryAry.h and ptrAry.c

- Dependencies:
  - ulCp

# ptrAry

Has the structure and functions to manage an array of
  c-strings.

The core of ptrAry is the str\_ptrAry structure, which
  holds the array of c-string pionters. After makeing
  your structure you can add strings, sort strings, or
  search for strings.

# Structure

## structure variables

The str\_ptrAry holds the array of c-strings and the
  supporting variables to manage you c-string array.

- Variables:
  1. strAry: arrray of c-string pionters
             (singed char \*\*) to have input strings
  2. lenAryUI: array with the length of every string in
               strAry (var 1)
  3. sizeAryUI: array with the maximum bytes each c-string
                in strAry (var 1) can get before needing
                to be resized
  4. sizeSL: maximum strings in strAry (var 1) before
     needing to resize
  5. lenSL: number of strings currently in strAry (var 1)

## structure setup functions

You can setup/initialize/free a str\_ptrAry structure
  using the blank, init, mk, resize, and free functions.

### quick structure make and setup

The quick way to get a str\_ptrAry is to using
  the `mk_str_ptrAry()`. This function takes in the number
  of c-strings you want to start out with. This function
  will create the structure (on heap), initialize, and
  then allocate memory for the structure.

The output is a pointer to a heap allocated str\_ptrAry
  structure or 0 for memory errors.

```
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);

/*code goes here*/
freeHeap_str_ptrAry(strAryHeapST); /*will go over later*/
strAryHeapST = 0;
```

### manually setting up structure

You can skip the `mk_str_ptrAry()` step, but in return
  you have to call the functions it calls behind the
  scenes.

First, declare your structure, then initialize your
  structure, finally allocate memory for your c-string
  array.

You can initialize your structure using
  the `init_str_ptrAry()` function. This sets all ponters
  to 0/null (so only call once). The input is a pointer
  to your str\_ptrAry structure to initialize.

Next, you can allocate (or reallocate) memory for the
  arrays using the `resize_str_ptrAry()` function.

- Input:
  1. str\_ptrAry structure pointer to allocate/reallocate
     memory
  2. number of c-strings for the structure to hold
- Output:
  - Returns: 0 for no errors
  - Returns: def\_memErr\_ptrAry for memory errors

```
struct str_ptrAry cStringAryStackST;
init_str_ptrAry(&cStringAryStackST);

if( resize_str_ptrAry(&cStringAryStackST, 10) )
   /*deal with memory errors*/

/*your code goes here*/

freeStack_str_ptrAry(&cStringAryStackST); /*next section*/
```

### freeing your structure

When finished with your str\_ptrAry structure, you need to
  free it. There are two functions to manage freeing.

Both functions take as input a pointer to your str\_ptrAry
  function. There is no output (beyond freeing).

The `freeHeap_str_ptrAry()` will free a heap allocated
  str\_ptrAry structure, such as returned from
  `mk_str_ptrAry()`. After calling `freeHeap_str_ptrAry`,
  make sure to set the pointer to 0.

```
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);

/*code goes here*/
freeHeap_str_ptrAry(strAryHeapST);
strAryHeapST = 0;
```

The `freeStack_str_ptrAry()` function will free the
  variables inside a str\_ptrAry structure and set it to
  an initialized state. Use this to free stack allocated
  structures.

```
struct str_ptrAry cStringAryStackST;
init_str_ptrAry(&cStringAryStackST);

if( resize_str_ptrAry(&cStringAryStackST, 10) )
   /*deal with memory errors*/

/*your code goes here*/

freeStack_str_ptrAry(&cStringAryStackST);
```

# Functions

These are functions you would use to blank, add strings,
  sort strings, or look up strings for a str\_ptrAry
  structure.

## blank structures

You can blank all c-strings (set first elemnt to null)
  with the `blank_str_ptrAry()` function. This will also
  set each lenAryUI element to 0.

- Input: str\_ptrAry pointer to structure to blank

```
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);

/*code goes here*/

blank_str_ptrAry(strAryHeapST); /*set all c-strings to 0*/

/*code after blank goes here*/

freeHeap_str_ptrAry(strAryHeapST);
strAryHeapST = 0;
```

## adding strings

You can add strings to a str\_ptrAry structure using
  the `add_str_ptrAry()` function. This function handels
  resizing automatically.

- Input:
  1. c-string to add (copy) to the c-string array
  2. str\_ptrAry structure pointer with c-string array
  3. index to add new c-string (input 1) at
- Output:
  1. Modifies: strAry to have new c-string (input 1) at
               index (input 3)
  2. Modifies: lenAryUI to have c-string (input 1) length
  3. Returns:
     - 0 for no errors
     - def\_memErr\_ptrAyr

```
signed char errSC = 0;
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);

errSC =
   add_str_ptrAry(
      (signed char *) "hellow world!",
      strAryHeapST,
      0 /*at index 0; so strAryHeapST->strAry[0]*/
   );
if(errSC)
   /*deal with memory error*/

errSC =
   add_str_ptrAry(
      (signed char *) "by world!",
      strAryHeapST,
      1 /*at index 1; so strAryHeapST->strAry[1]*/
   );
if(errSC)
   /*deal with memory error*/

errSC =
   add_str_ptrAry(
      (signed char *) "no world",
      strAryHeapST,
      0 /*at index 0; so overwrites "hellow world!"*/
   );
if(errSC)
   /*deal with memory error*/

freeHeap_str_ptrAry(strAryHeapST);
strAryHeapST = 0;
```

## sorting

You can sort your str\_ptrAry with the sort function
  (uses shell sort) or by adding strings in a sorted
  order (using a binrary search, but still need to move
  strings around, so O(n^2)).

### sorting during addition

To add strings in a sorted order you must use
  the `addSort_str_ptrAry()` function. You can not use
  `add_str_ptrAry()`, unless you sort the strings before
  calling `addSort_str_ptrAry()`. This will also mess up
  your index.

The one warning is this is a O(n^2) time do to shifting
  strings around.

- Input:
  1. c-string to add in
  2. str\_ptrAry structure (as pointer) to add c-string to
- Output:
  1. inserts in sorted position c-string and adjusts
     lenghts
  2. Returns:
     - index of added string
     - -1 for memory errors

```
signed char errSC = 0;
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);

errSC =
   addSort_str_ptrAry(
      (signed char *) "hellow world!",
      strAryHeapST
   );
if(errSC)
   /*deal with memory error*/

errSC =
   addSort_str_ptrAry(
      (signed char *) "by world!",
      strAryHeapST
   );
if(errSC)
   /*deal with memory error*/


freeHeap_str_ptrAry(strAryHeapST);
strAryHeapST = 0;
```

### sorting when finished

The other option is to add strings and then sort them
  later with the `sort_str_ptrAry()` function. The input
  is the structure to sort.

You can also keep an unsigned int array in sync with the
  strings using the `sortSync_str_ptrAry()` function.

- Input `sortSync_str_ptrAry()`:
  1. str\_ptrAry (as pointer) structure to sort
  2. unsigned int array to keep in sync


```
signed char errSC = 0;
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);
unsigned int *indexArrayUI[2] = {0, 1};

errSC =
   add_str_ptrAry(
      (signed char *) "hellow world!",
      strAryHeapST,
      0
   );
if(errSC)
   /*deal with memory error*/

errSC =
   add_str_ptrAry(
      (signed char *) "by world!",
      strAryHeapST,
      1
   );
if(errSC)
   /*deal with memory error*/

sortSync_str_ptrAry(strAryHeapST, indexArrayUI);
   /*this keeps the index array in sync
   `   index 0:
   `     - strAryHeapST "by world!
   `     - indexArrayUI 1
   `   index 1:
   `     - strAryHeapST "hellow world!
   `     - indexArrayUI 0
   */

/*if sorting strings, but not keeping array in sync*/
sort_str_ptrAry(strAryHeapST);

freeHeap_str_ptrAry(strAryHeapST);
strAryHeapST = 0;
```

## finding c-strings

You can find c-strings using the `find_str_ptrAry()`
  and the `findNoSort_str_ptrAry()` functions. Both
  functions take the same input and return the same
  output.

The difference is that `find_str_ptrAry()` does a binary
  search and needs a sorted array,
  while `findNoSort_str_ptrAry()` starts from the start
  and moves till the string is found or the end of the
  array.

- Input (both):
   1. str\_ptrAry structure (as pointer) to search for
      c-string int
   2. c-string to search for
- Output:
   - Returns:
     - index of found c-string
     - -1 if could not find the c-string

```
signed char errSC = 0;
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);
signed int indexSI = 0;

errSC =
   add_str_ptrAry(
      (signed char *) "hellow world!",
      strAryHeapST,
      0
   );
if(errSC)
   /*deal with memory error*/

errSC =
   add_str_ptrAry(
      (signed char *) "by world!",
      strAryHeapST,
      1
   );
if(errSC)
   /*deal with memory error*/

sort_str_ptrAry(strAryHeapST); /*sort the arrays*/

indexSI =
   find_str_ptrAry(
      strAryHeapST,
      (signed char *) "hellow world!"
   ); /*returns 1 (sort puts at 1)*/

indexSI =
   find_str_ptrAry(
      strAryHeapST,
      (signed char *) "bellow world!"
   ); /*returns -1 (not in array)*/

freeHeap_str_ptrAry(strAryHeapST);
strAryHeapST = 0;
```

## misc functions

### string comparison

You can compare two strings in a str\_ptrAry structure
  using the `cmp_str_ptrAry()` function.

- Input:
  1. str\_ptrAry structure (as pointer) with c-strings
  2. index of query (firsts string) to compare
  3. index of reference (second string) to compare
- Output:
  - Returns:
    - 0 if equal
    - > 0 if query (input 2) is greater
    - < 0 if reference (input 3) is greater

```
signed char errSC = 0;
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);
signed long resultSL = 0;

errSC =
   add_str_ptrAry(
      (signed char *) "hellow world!",
      strAryHeapST,
      0
   );
if(errSC)
   /*deal with memory error*/

errSC =
   add_str_ptrAry(
      (signed char *) "by world!",
      strAryHeapST,
      1
   );
if(errSC)
   /*deal with memory error*/

resultSL = cmp_str_ptrAry(strAryHeapST, 0, 1);
  /*returns > 0*/ 

resultSL = cmp_str_ptrAry(strAryHeapST, 0, 0);
  /*returns 0*/ 

freeHeap_str_ptrAry(strAryHeapST);
strAryHeapST = 0;
```

### string swapping

You can swap two strings in a str\_ptrAry structure
  using the `swap_str_ptrAry()` function.

- Input:
  1. str\_ptrAry structure (as pointer) with c-strings
  2. index of first string to swap
  3. index of second string to swap

```
signed char errSC = 0;
struct str_ptrAry *strAryHeapST = mk_str_ptrAry(10);
signed long resultSL = 0;

errSC =
   add_str_ptrAry(
      (signed char *) "hellow world!",
      strAryHeapST,
      0
   );
if(errSC)
   /*deal with memory error*/

errSC =
   add_str_ptrAry(
      (signed char *) "by world!",
      strAryHeapST,
      1
   );
if(errSC)
   /*deal with memory error*/

swap_str_ptrAry(strAryHeapST, 0, 1);
  /*moves "hellow world!" to index 1 and
  `  "by world!" to index 0
  */ 

freeHeap_str_ptrAry(strAryHeapST);
strAryHeapST = 0;
```

### insert index

You can find the index to insert a string in so that the
  array remains sorted using the `findClose_str_ptrAry()`
  function.

You must make sure your array is sorted before calling
  this function. It uses a binrary search.

- Input:
  1. str\_ptrAry structure (as pointer) with c-strings
  2. c-string to find insert location at
- Output:
  - Returns: best location to insert c-string (input 2)
             so that the array remains sorted

I would reccomend using `addSort_str_ptrAry()` instead,
  since that calls this function.
