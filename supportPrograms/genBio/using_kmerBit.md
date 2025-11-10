# Goal

Describe how to use the kmerBit library to convert
  nucleotide sequences into numeric kmers.

Files: kmerBit.h

Depenedencies: none

You should be familar with the ntTo2Bit library (see
  using\_ntTo2Bit.md). It is also recomended to be familar
  with the revNtTo2Bit library (see using\_revNtTo2Bit.md)
  and alnSet (see ../alnSet/using\_alnSet.md) libraries.

# Overview

The kmerBit.h file has one look up table and three macros.

The flow is to first build a mask for clearing old
  nucleotides from the kmer with `mkMask\_kmerBit()`. Then
  convert the target nucleotide using the look up tables
  in ntTo2Bit.h or revNtTo2Bit.h libraries. After that
  you can use `ntBitToKmer\_kmerBit()` function to add the
  nucleotide to the kmer.

As a rule you should look for anonymous bases and errors
  from the `ntTo2Bit.h` lookup tables. These functions
  assume that the input is a regular (`ATCG`) nucleotide.

If you have converted you sequences
  using the `seqToIndex\_alnSet()` function
  in `../genAln/alnSet.c/h`, then you can use
  the `ntIndexToKmer\_kmerBit()` function to add the
  nucleotide indexed to the kmer.

The other opton would be to use the lookup table to
  convert the indexed base to a ntTo2Bit.h code. Again,
  watch out for anonymous bases or errors.

# make a mask

The mask in the kmer bit function is designed to clear
  the old nucloetides out of the kmer. You can build it
  with the `mkMask\_kmerBit()` function. The input is the
  length of the kmer (max 32 nucleotides). The output is
  an unsigned long with the mask.


for example, to make a mask for 3mers, you would
  use `unsigned long maskUL = mkMask\_kmerBit(3);`.

# add nucleotide to kmer

Once you have the mask made you can then convert your
  nucleotide to a two bit value and then add it to your
  kmer with the `ntBitToKmer\_kmerBit();` function.

- Input
  1. nucleotdie in two bit (0 to 3) format to add
  2. kmer currently working on (unsigned long)
  3. the mask returned from `mkMask_kmerBit()`
- Output:
  - returns input 2 with the old nucleotide removed and
    the new nucleotide added int

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "ntTo2Bit.h"
#include "kmerBit.h"

int
main(
){
   signed char *seqStr =
      (signed char *) "ATCATGACTGACGATGCGATGC";
   signed char *baseStr = seqStr;
   signed int kmerOnSI = 0;

   /*since the kmer is less then 16 I am using an int*/
   unsigned int maskUI = mkMask_kmerBit(7);
   unsigned int kmerUI = 0;
   unsigned char ntUC = 0;
   signed int lenKmerSI = 0;

   while(*baseStr)
   { /*Loop: convert sequence to kmers*/
      ++kmerOnSI;
      ntUC = ntTo2Bit[(unsigned char) *baseStr];

      if(ntUC > 3)
      { /*If: had an anoymous base or invalid character*/
         kmerUI = 0;
         lenKmerSI = 0;
      } /*If: had an anoymous base or invalid character*/

      else
      { /*Else: have a valid base*/
         kmerUI = ntBitToKmer_kmerBit(ntUC,kmerUI,maskUI);

         if(lenKmerSI < 7)
            ++lenKmerSI;
         else
            printf("kmer %i is %u\n", kmerOnSI, kmerUI);
      } /*Else: have a valid base*/

      ++baseStr;
   } /*Loop: convert sequence to kmers*/

   return 0;
} /*main*/
```

# working with alnSet nucleotide index

If you are doing alignments using libraries from genAln,
  then it is likely you first converted you sequence into
  indexes with the `seqToIndex_alnSet()` function. In
  this case you have two options, if you know your
  sequence as no anonymous bases you can convert it using
  `ntIndexToKmer_kmerBit()`. However, if there is a chance
  it will, then you should use the lookup table and then
  call `ntBitToKmer_kmerBit()`.

## ntIndexToKmer\_kmerBit

This function converts and adds a nucleotide index from
  the `seqToIndex_alnSet()` function to a kmer. The one
  thing to be aware of is that this does not check for
  anonymous or invalid nucleotides.

- Input
  1. nucleotide index from `seqToIndex_alnSet()`
  2. kmer to add nucleotide to
  3. mask to remove old nucleotides with
- Output:
  - returns input 2 with the old nucleotide removed and
    the new nucleotide added int

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "../genAln/alnSet.h"
#include "kmerBit.h"

int
main(
){
   signed char seqStr[] =
      (signed char *) "ATCATGACTGACGATGCGATGC";
   signed char *baseStr = seqStr;
   signed int kmerOnSI = 0;

   /*since the kmer is less then 16 I am using an int*/
   unsigned int maskUI = mkMask_kmerBit(7);
   unsigned int kmerUI = 0;

   seqToIndex_alnSet(seqStr);

   while(*baseStr)
   { /*Loop: convert sequence to kmers*/
      ++kmerOnSI;
      kmerUI = ntIndexToKmer_kmerBit(ntUC,kmerUI,maskUI);

      if(lenKmerSI < 7)
         ++lenKmerSI;
      else
         printf("kmer %i is %u\n", kmerOnSI, kmerUI);

      ++baseStr;
   } /*Loop: convert sequence to kmers*/

   return 0;
} /*main*/
```

## index lookup table

The `alnNtTo_kmeBit` table converts a index nucleotide
  from the `seqToIndex_alnSet()` function to a 2 bit
  nucleotide that can be added to a kmer.

For anonymous bases it will return 4, for an invalid
  nucleotide it will return 8. So, make sure you check
  the return value before using.

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "../genAln/alnSet.h"
#include "kmerBit.h"

int
main(
){
   signed char seqStr[] =
      (signed char *) "ATCATGACTGACGATGCGATGC";
   signed char *baseStr = seqStr;
   signed int kmerOnSI = 0;

   /*since the kmer is less then 16 I am using an int*/
   unsigned int maskUI = mkMask_kmerBit(7);
   unsigned int kmerUI = 0;
   unsigned char ntUC = 0;
   signed int lenKmerSI = 0;

   seqToIndex_alnSet(seqStr);

   while(*baseStr)
   { /*Loop: convert sequence to kmers*/
      ++kmerOnSI;
      ntUC = alnNtTo_kmerBit[(unsigned char) *baseStr];

      if(ntUC > 3)
      { /*If: had an anoymous base or invalid character*/
         kmerUI = 0;
         lenKmerSI = 0;
      } /*If: had an anoymous base or invalid character*/

      else
      { /*Else: have a valid base*/
         kmerUI = ntBitToKmer_kmerBit(ntUC,kmerUI,maskUI);

         if(lenKmerSI < 7)
            ++lenKmerSI;
         else
            printf("kmer %i is %u\n", kmerOnSI, kmerUI);
      } /*Else: have a valid base*/

      ++baseStr;
   } /*Loop: convert sequence to kmers*/

   return 0;
} /*main*/
```
