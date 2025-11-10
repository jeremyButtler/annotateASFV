# Goal:

Describe how to use the kmerFun library to convert
  nucleotide sequences to kmer sequences.

Files: kmerFun.c and kmerFun.h

- Dependencies:
  - kmerBit.h

# Overview

The kmerFun library ony has one function
  (`seqToKmer_kmerFun()`) , which converts a nucleotide
  sequence to a kmer (at most 15mer) sequence.

There are two defined variabes. `def_endKmers_kmerFun`,
  which marks the end of the returned kmer array
  and `def_noKmer_kmerFun` which marks that a kmer had
  anonymous bases and so, no kmer was made.

The `seqToKmer_kmerFun()` function:

- Input:
  1. c-string with sequence to convert
  2. lenght of sequence to convert
  3. pointer to a singed int array (double pointer) to
     assign kmers to
     - resize/assign memory as needed
  4. signed int pointer with size of the kmer array
     (input 3)
     - changed when array is resized or assigned memory
  5. length of one kmer
- Output:
  - Returns number of kmers in the sequence or 0 if had
    memory errors
  - fills the kmer array (input 3) with kmers

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "kmerFun.h"

int
main(
){
   signed char *seqenceStr = "ATGCGATGCAGGTAAAAAGCTTTTT";
   signed int sequenceLenSI = 25; /*25 bases long*/

   unsigned char kmerLengthUC = 7;

   signed int *kmerHeapArySI = 0;
   signed int arraySizeSI = 0;
   signed int numberKmersSI = 0;
   signed int siKmer = 0;

   numberKmersSI =
      seqToKmer_kmerFun(
         sequenceStr,    /*sequence*/
         sequenceLenSI,  /*length of sequence*/
         &kmerHeapArySI, /*function will assign memory*/
         &arraySizeSI,   /*gets kmer array size*/
         kmerLengthUC
      );

   if(! numberKmersSI)
   { /*If: memory error*/
      fprintf(stderr, "ran out of memory\n");
      return 1;
   } /*If: memory error*/

   printf(
      "Sequence had %i kmers:", 
      numberKmersSI
   );

   /*print out the kmers in the sequence*/
   for(siKmer = 0; siKmer < numberKmersSI; ++siKmer)
      printf(" %i", kmerHeapArySI[siKmer]);

   printf("\n");

   free(kmerHeapArySI);
   kmerHeapArySI = 0;
   return 0;
} /*main*/
```
