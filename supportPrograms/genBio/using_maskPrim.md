# Goal

Describe how to use the maskPrim library to convert
  mask known primers coordinates in sam file entries
  (reads).

Files: maskPrim.h maskPrim.c

- Depenedencies:
  - samEntry
  - fileFun from genLib
  - shellSort from genLib
  - base10Str from genLib

- Hidden dependencies:
  - numToStr from genLib
  - ulCp from genLib
  - endLine.h from genLib
  - ntTo5Bit.h

The maskPrim library using the samEntry structure, so you
  should be familar with the samEntry library (see
  using\_samEntry.md).

# Overview

The maskPrim system first reads in coordinates from a tsv
  file. It then uses these coordinates to mask reads in
  a sam file. The coordinates are treated as if they were
  primer locations, since it is typically primer locations
  that are masked.

## How the system works:

As a rule you do not need to think about this. There are
  macros which are designed to handel this. Also, there
  is a read in function that will set this up
  automatically.

Primers are stored in an unsigned int array, were each
  primer takes up one unsigned int.

1. First bit encodes the direction
   - 1 for reverse complement primer
   - 0 for foward primer
2. Second bit encodes if the primer should always be
   treated as a pair
   - 1 always paired; make sure both
   - 0 treat each primer as unpaired (independent) primer
3. The remaing bits store the index of the primers pair if
   the paired bit is set to 1

For example; say primers `0x00000001`, `0x00001110`, and
  `0x00001011` were in an array.

```
0x00000001 primer is reverse complemnt, but unpaired
         ^-1 here means reverse complement

        v----means primer is paired
0x00001110 primer is forward and paired to primer 3
      ^^
      |+______ index of paired primer: 0x000011 = index 3

        +----means primer is paired
        Vv----means primer is reverse complement
0x00001011 primer is reverse complement and is paired to
      ^^     primer number 2
      |+______ index of paired primer: 0x000010 = index 2
```

The primer array is paired up with a primer start array
  and a primer end array. These arrays hold the primer
  starting and ending coordinates.

So you have three arrays, the primer data array
  (direction, if paired, index of pair), the starting
  coordinate array of each primer, and the ending
  coordinate array of each primer. These arrays always
  need to be kept in synce, so if you sort one or swap
  elements in one you need to sort/swap elements in all
  three.

As a rule, you should not have to worry about messing with
  primers. This is here so you can understand how the data
  is stored.

## reading in coordinates

You can read in a tsv (tab or space deliminated) file with
  primer coordinates with the `getCoords_maskPrim()`
  function. The tsv file format is:

First row is the header and the remaining rows are for the
  primers.

- Column 1: is ignored
  - suggest primer name or reference id
- Column 2: is ignored
  - suggest primer name or reference id
- Column 3: if the primers are paried or separate
  - you can use `True` and `False`, `1` and `0`, or
    `+` and `-`
- Column 4: starting coordinate of forward primer
  - use `NA` for no forward primer
- Column 5: ending coordinate of forward primer
  - use `NA` for no forward primer
- Column 6: starting coordinate of reverse primer
  - use `NA` for no reverse primer
- Column 7: ending coordinate of reverse primer
  - use `NA` for no reverse primer

Example with three masking coordinates. The frist one is a
  single 20 base position to maske. The second requires
  both positions to mask. The third is a single reverse
  position (no idea why you would do this).

```
id reference paired forStart forEnd revStart revEnd
1  NC000962  False  10       30     NA       NA
2  NC000962  True   100      122    500      521
3  NC000962  False  NA       NA     1000     1018
```

- Input:
  1. c-string with file name to open
  2. unsigned int array pointer
     (`unsigned int **startAryUI=0;`) to hold starting
     coordinates of each primer
  3. unsigned int array pointer
     (`unsigned int **endAryUI=0;`) to hold ending
     coordinates of each primer
  4. unsigned int array to hold the direction, paired
     status, and index of paired primer
     (`unsigned int **primersAryUI`)
  5. signed long pointer to get errors
     - set to 0 for no errors
     - set to def\_memErr\_maskPrim for memory errors
     - set to `lineX << 8 | def_fileErr_maskPrim` for
       file error on line number x
     - set to `def_emptyFileErr_maskPrim` if file has
       nothing
- Output:
  - Modifes: inputs 2, 3, 4, and 5 (how is  mentioned in
    input section)
    - the arrays are always sorted by starting coordinate
  - Returns: number of primers extracted
  - Returns: 0 for an error

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "maskPrim.h"
#include "samEntry.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char primersStr = (signed char *) "primers.tsv";

   signed int numberPrimersSI = 0;
   unsigned int *primersAryUI = 0;
   unsigned int *startAryUI = 0;
   unsigned int *endAryUI = 0;

   signed long errorSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: get coordinates
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   numberPrimersSI =
      getCoords_maskPrim(
         primersStr,
         &startAryUI,
         &endAryUI,
         &primersAryUI,
         &errorSL
      );

   if(! numberPrimersSI)
   { /*If: had an error*/
      if(errorSL == def_memErr_maskPrim)
         goto memErr_main_sec04;

      else if(errorSL == def_emptyFileErr_maskPrim)
      { /*Else If: nothing in file*/
         fprintf(stderr, "%s is empty\n", primersStr);
         goto fileErr_main_sec04;
      } /*Else If: nothing in file*/

      else if(errorSL & def_fileErr_maskPrim)
      { /*Else If: nothing in file*/
         fprintf(
            stderr,
            "error on line %li in %s\n",
            errorSL >> 8,
            primersStr
         );

         goto fileErr_main_sec04;
      } /*Else If: nothing in file*/
   } /*If: had an error*/

   printf("read in primers from %s\n", primersStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSL = 0;
   goto ret_main_sec04;

   fileErr_main_sec04:;
      errorSL = 2;
      goto ret_main_sec04;

   memErr_main_sec04:;
      fprintf(stderr, "ran out of memory\n");
      errorSL = 1;
      goto ret_main_sec04;

   ret_main_sec04:;
      if(primersAryUI)
         free(primersAryUI);
      primersAryUI = 0;

      if(startAryUI)
         free(startAryUI);
      startAryUI = 0;

      if(endAryUI)
         free(endAryUI);
      endAryUI = 0;

      return errorSL;
} /*main*/
```

## mask positions

Once you have the primer coordinates you can then mask
  the locations in the reads using the `maskPrim()`
  function.

- Input:
  1. samEntry structure with read to mask
  2. unsigned int array with masking start coordinates
  3. unsigned int array with masking end coordinates
  4. unsigned int array with masking primers (flags)
  5. number of primers/positions to mask
  6. number of bases in to look for pimers
     - -1 is no limit, mask any position in read
     - 0 is only if first/last base in read
     - x > 0 is x base from start or before end
  7. character to mask with (suggest `N`)
- Modifies:
  1. input 1 to have coordinates masked
- Returns:
  1. 0 for success
  2. -1 if no masking sites in read
  3. (index << 1) | 1 if a masking site was present, but
     was out of range (input 6 sets range)

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "maskPrim.h"
#include "samEntry.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char primersStr = (signed char *) "primers.tsv";

   signed int numberPrimersSI = 0;
   unsigned int *primersAryUI = 0;
   unsigned int *startAryUI = 0;
   unsigned int *endAryUI = 0;

   signed long errorSL = 0;

   struct samEntry samStackST;
   FILE *samFILE = fopen("reads.sam", "r");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: get coordinates
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! samFILE)
   { /*If: no sam file input*/
      fprintf(stderr, "reads.sam does not exist\n");
      goto fileErr_main_sec04;
   } /*If: no sam file input*/

   init_samEntry(&samStackST);
   if( setup_samEntry(&samStackST) )
      goto memErr_main_sec04;

   numberPrimersSI =
      getCoords_maskPrim(
         primersStr,
         &startAryUI,
         &endAryUI,
         &primersAryUI,
         &errorSL
      );

   if(! numberPrimersSI)
   { /*If: had an error*/
      if(errorSL == def_memErr_maskPrim)
         goto memErr_main_sec04;

      else if(errorSL == def_emptyFileErr_maskPrim)
      { /*Else If: nothing in file*/
         fprintf(stderr, "%s is empty\n", primersStr);
         goto fileErr_main_sec04;
      } /*Else If: nothing in file*/

      else if(errorSL & def_fileErr_maskPrim)
      { /*Else If: nothing in file*/
         fprintf(
            stderr,
            "error on line %li in %s\n",
            errorSL >> 8,
            primersStr
         );

         goto fileErr_main_sec04;
      } /*Else If: nothing in file*/
   } /*If: had an error*/

   printf("read in primers from %s\n", primersStr);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: mask reads
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSL = get_samEntry(&samStackST, samFILE);

   if(! errorSL)
      ;
   else if(errorSL == def_memErr_samEntry)
      goto memErr_main_sec04;
   else
   { /*Else: file error*/
      fprintf(stderr, "file error reading reads.sam\n");
      goto fileErr_main_sec04;
   } /*Else: file error*/


   while(! errorSL)
   { /*Loop: mask reads*/
      if(! samStackST.seqStr[0])
         goto nextRead_main_sec03;
      if(samStackST.seqStr[0] == '*')
         goto nextRead_main_sec03;

      errorSL =
         maskPrim(
            &samStackST,  /*has read to mask*/
            startAryUI,   /*start coordinates of primers*/
            endAryUI,     /*end coordinates of primers*/
            primersAryUI, /*primer information*/
            numberPrimersSI, /*number primers user input*/
            -1,           /*mask at any position*/
            'N'           /*mask bases to N*/
         );

      p_samEntry(&samStackST, 0, stdout);

      nextRead_main_sec03:;
         errorSL = get_samEntry(&samStackST, samFILE);
   } /*Loop: mask reads*/


   if(errorSL != def_EOF_samEntry)
      ;
   else if(errorSL == def_memErr_samEntry)
      goto memErr_main_sec04;
   else
   { /*Else: file error*/
      fprintf(stderr, "file error reading reads.sam\n");
      goto fileErr_main_sec04;
   } /*Else: file error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSL = 0;
   goto ret_main_sec04;

   fileErr_main_sec04:;
      errorSL = 2;
      goto ret_main_sec04;

   memErr_main_sec04:;
      fprintf(stderr, "ran out of memory\n");
      errorSL = 1;
      goto ret_main_sec04;

   ret_main_sec04:;
      if(samFILE)
         fclose(samFILE);
      samFILE = 0;

      if(primersAryUI)
         free(primersAryUI);
      primersAryUI = 0;

      if(startAryUI)
         free(startAryUI);
      startAryUI = 0;

      if(endAryUI)
         free(endAryUI);
      endAryUI = 0;

      return errorSL;
} /*main*/
```

## other functions

There are three kinds of other functions, the first are to
  add primer direction, paired status, or index of paired
  primer to the primer array. The second are to get the
  direction, paired status, and paired primers index
  from the primer array. The thrid is a function to sort
  the primers.

### add data to  primer array

These functions are really macros, so be carefull.

For these functions, primersList is an unisgned int array
  of primers. indexUI is the index of one primer in the
  array. pairIndexUI is the index of a second primer the
  indexUI primer is paired to.

These are functions for setting the primer values in an
  array.

1. `blankFlag_maskPrim()` sets a single primers direction,
  paired status, and mate index to zero
   - input is unsigned int with primer to blank
   - `blankFlag_maskPrim(primerList[indexUI])`
2. `setDir_maskPrim()` sets a single primers direction to
  the input value
   - for reverese `setDir_maskPrim(primerList[indexUI], 1)`
   - for forward `setDir_maskPrim(primerList[indexUI], 0)`
3. `setPair_maskPrim()` sets the primer to have a paired
   or unpaired state
   - for paired `setPair_maskPrim(primerList[indexUI], 1)`
   - for unpaired `setPair_maskPrim(primerList[indexUI], 0)`
4. `setMateIndex_maskPrim()` sets the index for the primer
   a primer is paired to
   - `setMateIndex_maskPrim(primerList[indexUI], pairIndexUI)`

### get data from primer array

These are functions for getting the values for a primer

These functions are really macros, so be carefull.

For these functions, primersList is an unisgned int array
  of primers. indexUI is the index of one primer in the
  array.

1. `getDir_maskPrim()` gets the direction (as 1 for
   reverse or 0 for forward) from a primer
   - `getDir_maskPrim(primerList[indexUI])`
2. `getPair_maskPrim()` gets if the primer is paired (1)
   or unpaired (0)
   - `getPair_maskPrim(primerList[indexUI])`
3. `getMateIndex_maskPrim()` gets the index of the primer
   a primer is paired with
   - `getMateIndex_maskPrim(primerList[indexUI])`

### sorting primers

Finally there is a sorting function
  named `sortCoords_maskPrim()`. This is called at the end
  of `getCoords_maskPrim()`. So, only call this if you are
  using something else.

- Input:
  1. unsigned int array with primer starting coordinates
  2. unsigned int array with primer ending coordinates
  3. unsigned int array with primer data (direction,
     paired, and index of paired primer)
  4. number of primers in the array
- Modifies:
  1. input arrays to be ordered by starting coordinates
     array (least to greates sort)

```
   numberPrimersSI =
      getCoords_maskPrim(
         primersStr,
         &startAryUI,
         &endAryUI,
         &primersAryUI,
         &errorSL
      );

   /*did something to mess up sorting order*/

   sortCoords_maskPrim(
      startAryUI,
      endAryUI,
      primersAryUI,
      numberPrimersSI
   );
```
