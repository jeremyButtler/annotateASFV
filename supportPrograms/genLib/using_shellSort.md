# Goal:

Describing how to use the shellShort library to sort
  arrays and search sorted arrays.

Files: shellShort.h and shellShort.c

# shellSort

The shellSort library stores functions to sort each
  numeric datatype with shell sort. It also includes two
  binrary search functions.

Because there is a set of three functions for each
  datatype (except doubles) there are a lot of duplicate
  functions. So, I will list the general example and
  unique functions.

## general functions

- `<data_type>Search_shellSort`: finds an item in sorted
  array with a binary search
  - Input:
    1. sorted array to search for input
    2. number to search for
    3. length of the array (input 1) (index 1)
  - Output:
    - Returns:
      - index of item
      - -1 if could not find the item
- `<data_type>Range_shellSort`: finds the first item
  in a sorted array that falls within the input range with
  a binary search
  - Input:
    1. sorted array to search for input
    2. start of range
    3. end of range
    4. length of the array (input 1) (index 1)
  - Output:
    - Returns:
      - index of fist item in range
      - -1 if could not find the item
- `<data_type>_shellSort`: sorts numeric array using
  shell sort
  - Input:
    1. array to sort
    2. first index in array to sort (use 0)
    3. last index in array to sort (index 0;
       use length - 1)

| data type      | short hand (in function name) |
|:---------------|:------------------------------|
| unsigned long  | ul or UL                      |
| unsigned int   | ui or UI                      |
| unsigned short | us or US                      |
| unsigned char  | uc or UC                      |
| signed long    | sl or SL                      |
| signed int     | si or SI                      |
| signed short   | ss or SS                      |

Table: short hand for datatypes in shellSort.h/c functions

**Example**

I am using the unsigned long functions as an example.

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "shellSort.h"

int
main(
){
   unsigned long ulAry[10] =
      {1, 3, 5, 2, 4, 7, 9, 0, 6, 8}; /*is 0-9*/

   ul_shellSort(ulAry, 0, 9); /*sort array*/

   printf(
      "index of 6: %ul\n",
      ulSearch_shellSort(ulAry, 6, 9)
   ); /*prints 6*/

   printf(
      "range of 3 to 6: %ul\n",
      ulRange_shellSort(ulAry, 3, 6, 9)
   ); /*prints 3 (fist item between 3 and 6)*/
} /*main*/
```

## specialized functions

These are functions I use for unique events. They are one
  off and so are only for one data type.

- `dbl_shellSort`: sorts an array of doubles with shell
  sort (same input as `<data_type>_shellShort`)
  - Input:
    1. double array to sort
    2. first index in array to sort (use 0)
    3. last index in array to sort (index 0;
       use length - 1)
- `uiTwinSort_shellSort`: sorts the first unsigned int
  array and keeps a second unsigned int array in sync
  with the first
  - Input:
    1. unsigned int array to sort
    2. 2nd unsigned int array to keep in sync with first
    3. first index in array to sort (use 0)
    4. last index in array to sort (index 0;
       use length - 1)

**Example (twin sort)**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "shellSort.h"

int
main(
){
   signed int indexSI = 0;
   signed int oldIndexSI = 0;

   unsigned int uiAry[10] =
      {1, 3, 5, 2, 4, 7, 9, 0, 6, 8}; /*is 0-9*/
   unsigned int indexAryUI[10] =
      {0, 1, 2, 3, 4, 5, 6, 7, 8, 9};
      /*index of each item*/

   uiTwinSort_shellSort(uiAry, indexAryUI, 0, 9);
     /*sort array; but keep index array in sync with
     `  uiAry (number in indexAryUI stays with its number
     `  in uiAry)
     */

   /*I have the original index's preserved so I can
   `  unsort the arrays to the original order
   */
   for(indexSI = 0; indexSI < 10; ++indexSI)
   { /*Loop: unsort the arrays*/
      oldIndexSI = indexAryUI[indexSI];

      /*swapping elements around*/
      uiAry[indexSI] ^= uiAry[oldIndexSI];
      uiAry[oldIndexUI] ^= uiAry[indexSI];
      uiAry[indexSI] ^= uiAry[oldIndexSI];

      /*swapping elements around*/
      indexAryUI[indexSI] ^= indexAryUI[oldIndexSI];
      indexAryUI[oldIndexUI] ^= indexAryUI[indexSI];
      indexAryUI[indexSI] ^= indexAryUI[oldIndexSI];
   } /*Loop: unsort the arrays*/
} /*main*/
```
