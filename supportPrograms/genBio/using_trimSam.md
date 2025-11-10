# Goal

Describe how to use the trimSam library to trim reads
  in a sam file.

Files: trimSam.h trimSam.c

- Depenedencies:
  - samEntry
  - ulCp from genLib

- Hidden dependencies:
  - fileFun from genLib
  - endLine.h from genLib
  - numToStr from genLib
  - base10Str from genLib
  - ntTo5Bit.h

The trimSam library uses the samEntry structure, so you
  should be familar with the samEntry library (see
  using\_samEntry.md).

# Overview

The trimSam system allows you to trim reads in a sam file.

- Currently you can:
  1. trim off softmasking at ends (`seq_trimSam()`)
  2. trim by coordinates (`coords_trimSam()`)
  3. trim x bases of the ends of each read
    (`ends_trimSam()`)
  4. also includes a function to trim softmasked bases
     from the ends of all reads in a
      file (`trimReads_trimSam()`)

## 1 trim soft masked bases

You can trim soft masked bases off the ends of a reads
  using the `seq_trimSam()` function. The input is a
  samEntry structure pointer.

- Output:
  1. 0 if trimmed of soft masked bases (input is modified)
  2. def\_header\_trimSam if input a sam file header
     instead of a read
  3. def\_noMap\_trimSam if read was unmapped
  4. def\_noSeq\_trimSam if read had no sequence

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "trimSam.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   struct samEntry readStackST;
   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   init_samEntry(&readStackST);
   if( setup_samEntry(&readStackST) )
      goto memErr_main_sec04;
   
   samFILE = fopen("reads.sam", "r");
   if(! samFILE)
      goto fileErr_main_sec04;
   
   errorSC = get_samEntry(&readStackST, samFILE);
   
   if(! errorSC)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec04;
   else
      goto fileErr_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: trim reads
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   while(! errorSC)
   { /*Loop: add reads to the consnesus*/
     
      /*trim reads*/
      errorSC = seq_trimSam(&readStackST);

      /*make sure it is worth printing out*/
      if(errorSC == def_noMap_trimSam)
         ; /*read was unmapped, why print?*/
      else if(errorSC == def_noSeq_trimSam)
         ; /*no sequence, why print?*/
      else
         p_samEntry(&samStackST, 0, stdout);
         /*could trim or was a header*/

      errorSC = get_samEntry(&readStackST, samFILE);
   } /*Loop: add reads to the consnesus*/

   if(errorSC != def_EOF_samEntry)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec04;
   else
      goto fileErr_main_sec04;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto ret_main_sec04;

   memErr_main_sec04:;
      errorSC = 1;
      fprintf(stderr, "memory error\n");
      goto ret_main_sec04;

   fileErr_main_sec04:;
      errorSC = 2;
      fprintf(stderr, "file error with reads.sam\n");
      goto ret_main_sec04;

   ret_main_sec04:;
      if(! samFILE)
         ;
      else if(samFILE == stdin)
         ;
      else if(samFILE == stdout)
         ;
      else if(samFILE == stderr)
         ;
      else
         fclose(samFILE);
      samFILE = 0;
      
      freeStack_samEntry(&readStackST);
      return errorSC;
} /*main*/
```

## 2 trim by coordinates

You can trim reads by there mapping coordinates onto a
  reference using the `coords_trimSam()` function.

- Input:
  1. samEntry structure pointer with read to trim
  2. first reference base to trim at
  3. last reference base to trim at
- Output:
  1. 0 if trimmed of soft masked bases (input is modified)
  2. def\_header\_trimSam if input a sam file header
     instead of a read
  3. def\_noMap\_trimSam if read was unmapped
  4. def\_noSeq\_trimSam if read had no sequence
  5. def\_rangeErr\_trimSam if read had no bases between
     the start (input 2) and end (input 3)

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "trimSam.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   struct samEntry readStackST;
   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   init_samEntry(&readStackST);
   if( setup_samEntry(&readStackST) )
      goto memErr_main_sec04;
   
   samFILE = fopen("reads.sam", "r");
   if(! samFILE)
      goto fileErr_main_sec04;
   
   errorSC = get_samEntry(&readStackST, samFILE);
   
   if(! errorSC)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec04;
   else
      goto fileErr_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: trim reads
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   while(! errorSC)
   { /*Loop: add reads to the consnesus*/
     
      /*trim reads*/
      errorSC = coord_trimSam(&readStackST, 50, 1000);
         /*trim reads to be between reference position 50
         `  and 1000
         */

      /*make sure it is worth printing out*/
      if(errorSC == def_noMap_trimSam)
         ; /*read was unmapped, why print?*/
      else if(errorSC == def_noSeq_trimSam)
         ; /*no sequence, why print?*/
      else if(errorSC == def_rangeErr_trimSam)
         ; /*read has no bases between reference positions
           `  50 and 1000
           */
      else
         p_samEntry(&samStackST, 0, stdout);
         /*could trim or was a header*/

      errorSC = get_samEntry(&readStackST, samFILE);
   } /*Loop: add reads to the consnesus*/

   if(errorSC != def_EOF_samEntry)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec04;
   else
      goto fileErr_main_sec04;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto ret_main_sec04;

   memErr_main_sec04:;
      errorSC = 1;
      fprintf(stderr, "memory error\n");
      goto ret_main_sec04;

   fileErr_main_sec04:;
      errorSC = 2;
      fprintf(stderr, "file error with reads.sam\n");
      goto ret_main_sec04;

   ret_main_sec04:;
      if(! samFILE)
         ;
      else if(samFILE == stdin)
         ;
      else if(samFILE == stdout)
         ;
      else if(samFILE == stderr)
         ;
      else
         fclose(samFILE);
      samFILE = 0;
      
      freeStack_samEntry(&readStackST);
      return errorSC;
} /*main*/
```

## 3 trim by number

You can trim the last and first bases off reads with
  the `ends_trimSam()` function.

- Input:
  1. samEntry structure pointer with read to trim
  2. number of bases to trim from the start of the read
  3. number of bases to trim from the end of the read
  4. how strict to trim
     - 0: trim untill a match (`=` or `M`) is found
     - 1: apply 0 for the ending bases, but for the
          starting bases trim of the exact amount
     - 2: apply 0 for the starting bases, but for the
          ending bases trim of the exact amount
     - 3: trim of the exact number of input bases from
          the start and end
       - think of it as (1 | 2 or 1 combined with 2)
- Output:
  1. 0 if trimmed of soft masked bases (input is modified)
  2. def\_rangeErr\_trimSam if trimmed all bases from the
     the read
     - this likely includes headers, but not sure
  3. def\_noSeq\_trimSam if read had no sequence

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "trimSam.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   struct samEntry readStackST;
   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   init_samEntry(&readStackST);
   if( setup_samEntry(&readStackST) )
      goto memErr_main_sec04;
   
   samFILE = fopen("reads.sam", "r");
   if(! samFILE)
      goto fileErr_main_sec04;
   
   errorSC = get_samEntry(&readStackST, samFILE);
   
   if(! errorSC)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec04;
   else
      goto fileErr_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: trim reads
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   while(! errorSC)
   { /*Loop: add reads to the consnesus*/
      if(readStackST.extraStr[0] == '@')
      { /*If: this is a header*/
         p_samEntry(&readStackST, 0, stdout);
         goto nextRead_main_sec03;
      } /*If: this is a header*/
     
      /*trim reads*/
      errorSC = ends_trimSam(&readStackST, 10, 20, 0);
         /*trim 10 bases from start and 20 bases from end
         `  of the read. The 0 is trim untill match found
         */

      /*make sure it is worth printing out*/
      if(errorSC == def_noSeq_trimSam)
         ; /*no sequence, why print?*/
      else if(errorSC == def_rangeErr_trimSam)
         ; /*completely trimmed the read*/
      else
         p_samEntry(&samStackST, 0, stdout);
         /*could trim or was a header*/

      nextRead_main_sec03:;
         errorSC = get_samEntry(&readStackST, samFILE);
   } /*Loop: add reads to the consnesus*/

   if(errorSC != def_EOF_samEntry)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec04;
   else
      goto fileErr_main_sec04;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto ret_main_sec04;

   memErr_main_sec04:;
      errorSC = 1;
      fprintf(stderr, "memory error\n");
      goto ret_main_sec04;

   fileErr_main_sec04:;
      errorSC = 2;
      fprintf(stderr, "file error with reads.sam\n");
      goto ret_main_sec04;

   ret_main_sec04:;
      if(! samFILE)
         ;
      else if(samFILE == stdin)
         ;
      else if(samFILE == stdout)
         ;
      else if(samFILE == stderr)
         ;
      else
         fclose(samFILE);
      samFILE = 0;
      
      freeStack_samEntry(&readStackST);
      return errorSC;
} /*main*/
```

## 4 trim all reads in sam file

You can use the `trimReads_trimSam()` function to remove
  the starting and ending softmaksed bases from all reads
  in a sam file. Not very usefull.

- Input:
  1. FILE pionter to sam file with reads to trim
  2. FILE pionter to file to save trimmed reads to
  3. 1 or 0
     - 1: if want to keep unmapped reads
     - 0: if want to remove umapped reads
- Output:
  1. Prints trimmed to the output file (input 2)
  2. Returns:
     - 0 for success
     - def\_memErr\_samEntry for memory errors

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "trimSam.h"

int
main(
){
   signed char errorSC = 0;
   FILE *samFILE = 0;

   samFILE = fopen("reads.sam", "r");

   if(samFILE)
   { /*If: could open the file*/
      errorSC = trimReads_trimSam(samFILE, stdout, 0);
      /*0 means remove unmapped reads*/

      fclose(samFILE);
      samFILE = 0;

      if(errorSC)
         fprintf(stderr, "memory error\n");
   } /*If: could open the file*/

   else
   { /*Else: could not open the file*/
      fprintf(stderr, "could not open reads.sam\n");
      errorSC = 1;
   } /*Else: could not open the file*/

   return errorSC;
} /*main*/
```
