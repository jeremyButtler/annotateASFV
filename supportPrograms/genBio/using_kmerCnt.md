# Goal

Describe how to use the kmerCnt library to get the number
  of kmers in a sequence and find the number of matching
  kmers between two sequences.

Files: kmerCnt.c kmerCnt.h

- Dependencies:
  - seqST
  - genMath.h from genLib
  - ntTo2Bit.h
- Hidden Dependencies:
  - ulCp from genLib
  - fileFun from genLib
  - endLine from genLib

You should be familiar with the seqST library (seqST.c/h
  and using\_seqST.md).

# Overview:

The kmerCnt library is here to count the number of kmers
  that a sequence has. It can also find the number of
  matching kmers between two sequences. This is not
  designed to be memory efficent

# Structure:

The kmerCnt library uses the kmerCnt structure.

- Variables:
  - lenKmerUC: number of bases in a single kmer
    - max is 16 kmers
  - maxKmersUI: maximum number of kmers in a table
    (lenKmerUC * lenKmerUC)
  - forKmerArySI: signed int array with the forward
    sequences kmers
    - the size of the array is stored in forKmersUI
    - ends with def\_endKmers\_kmerCnt
    - kmers with anonymous bases are set to
      def\_noKmer\_kmerCnt
  - revKmerArySI: signed int array with the reverse
    complement sequences kmers
    - the size of the array is stored in revKmersUI
    - ends with def\_endKmers\_kmerCnt
    - kmers with anonymous bases are set to
      def\_noKmer\_kmerCnt
  - forSeqST: seqST structure pointer with the forward
    sequence
  - revSeqST: seqST structure pointer with the reverse
    complement sequence

# defined variables:

There are two types of defined variables, error return
  values and ones dealing with the kmer arrays in the
  kmerCnt structure.

- def\_noMatch\_kmerCnt:
  - no matches kmers found in both sequences
- def\_fileErr\_kmerCnt:
  - error reading in a fastx file
- def\_memErr\_kmerCnt:
  - memory error (most often ran out of memory)
- def\_noKmer\_kmerCnt:
  - the kmer had anonymous, invlaid character, or there
    was no kmer in the sequence
- def\_endKmers\_kmerCnt:
  - end of the kmer array

# Flow:

First create, initialize, and then setup a kmerCnt
  structure. Then add a seqeunce to the kmerCnt structure.
  After that, convert a target sequence to a kmer array
  and then compare the kmerCnt structure and target
  sequence number of matches. Finally, when finished free
  your kmerCnt structure and kmer array.

You can also convert a fasta file to an array of kmerCnt
  structures.

# Initialization and setup

You can initialize a kmerCnt structure using
  the `init_kmerCnt()` function. The input is a pointer to
  the kmerCnt structure to initialize. This function
  should only be called once.

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "kmerCnt.h"

int
main(
){
   struct kmerCnt kmerCountsStackST;
   init_kmerCnt(&kmerCountsStackST);

   freeStack_kmerCnt(&kmerCountsStackST);
} /*main*/
```

After initialization you can setup the memory in the
  kmerCnt structure using the `setup_kmerCnt()` function.
  Ths input is a pointer to the kmerCnt structure to
  setup and the length of one kmer. The ouput is 0 for
  success and def\_memErr\_kmerCnt for memory errors.

You should only call `setup_kmerCnt()` once per structer,
  or at least call `freeStack_kmerCnt()` to free the
  allocated memory first. The input is the kmerCnt
  structure to setup and the kmer length to use.

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "kmerCnt.h"

int
main(
){
   signed int errorSI = 0;
   struct kmerCnt kmerCountsStackST;
   init_kmerCnt(&kmerCountsStackST);

   if( setup_kmerCnt(&kmerCountsStackST, 3) )
      goto memoryError_main;


   errorSI = 0; /*no error*/
   goto cleanUp_main;

   memoryError_main:;
      fprintf(stderr, "ran out of memory\n");
      errorSI = 1;
      goto cleanUp_main;

   cleanUp_main:;
      freeStack_kmerCnt(&kmerCountsStackST);
      return errorSI;
} /*main*/
```

# converting sequences to kmers

## adding sequences to a kmerCnt structure

After initialization, you can add sequences to a kmerCnt
  structure using the `addSeq_kmerCnt()` function.

- Input:
  1. kmerCnt structure pointer to add the sequence to
  2. seqST structure pointer with the sequence to add
- Output:
  - 0 for no errors
  - def\_memErr\_kmerCnt for memory errors

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "seqST.h"
#include "kmerCnt.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec01: variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int errorSI = 0;
   struct seqST sequenceStackST;
   struct kmerCnt kmerCountsStackST;
   FILE *sequenceFILE = fopen("sequences.fasta", "r");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec02: initialize, setup, and get sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_kmerCnt(&kmerCountsStackST);
   init_seqST(&sequenceStackST);

   if( setup_kmerCnt(&kmerCountsStackST, 3) )
      goto memoryError_main_sec04;

   errorSI = getFa_seqST(sequenceFILE, &sequenceStackST);

   if(! errorSI)
      ; /*no errors*/
   else if(errorSI == def_memErr_seqST)
      goto memoryError_main_sec04;
   else
      goto fileError_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec04: return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSI = 0; /*no error*/
   goto cleanUp_main_sec04;

   memoryError_main_sec04:;
      fprintf(stderr, "ran out of memory\n");
      errorSI = 1;
      goto cleanUp_main_sec04;

   fileError_main_sec04:;
      fprintf(stderr, "error reading sequence.fasta\n");
      errorSI = 2;
      goto cleanUp_main_sec04;

   cleanUp_main_sec04:;
      freeStack_kmerCnt(&kmerCountsStackST);
      freeStack_seqST(&sequenceStackST);
      return errorSI;
} /*main*/
```

## adding a fasta file to kmerCnt structure

You can add all sequences in a fasta file to a kmerCnt
  structure array using the `faToKmerCnt_kmerCnt()`
  function. However, you will need to free the array
  with the array freeing functions.

You would do this if you planed to compare reads to a
  set of references.

- Input:
  1. c-string with the fasta to get seqeucnes from
  2. length of one kmer
  3. signed integer pointer to get the number of sequences
     in the fasta file
  4. signed character pointer to get errors
- Output:
  - 0 if had an error
    - input 4 gets set to
      - def\_fileErr\_kmerCnt for file errors
      - def\_memErr\_kmerCnt for file errors
  - heap allocated array of kmerCnt structures with kmer
    counts for each sequence

You can free the returened kmerCnt array using the
  `freeHeapAry_kmerCnt()` function.

- Input:
  1. the kmerCnt array to free
     - return value from `faToKmerCnt_kmerCnt()`
  2. the number of structures in the kmerCnt array
     - input 3 in `faToKmerCnt_kmerCnt()`
- Output:
  - nothing, but remember to set input 1 to 0/null

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "seqST.h"
#include "kmerCnt.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec01: variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   signed int numberSequencesSI = 0;
   signed int indexSI = 0;
   struct kmerCnt *kmerCountsAryHeapST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec02: initialize, setup, and get sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   kmerCountsAryHeapST = 
      faToKmerCnt_kmerCnt(
         (signed char *) "sequences.fasta",
         7, /*I want 7mers*/
         &numberSequencesSI,
         &errorSC
      );

   if(kmerCountsAryHeapST)
      ; /*no errors*/
   else if(errorSC == def_memErr_kmerCnt)
      goto memoryError_main_sec04;
   else
      goto fileError_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec03: do something with kmer counts
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(indexSI =0; indexSI < numberSequencesSI; ++indexSI)
      /*do something with your kmers*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec04: return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0; /*no error*/
   goto cleanUp_main_sec04;

   memoryError_main_sec04:;
      fprintf(stderr, "ran out of memory\n");
      errorSc = 1;
      goto cleanUp_main_sec04;

   fileError_main_sec04:;
      fprintf(stderr, "error reading sequence.fasta\n");
      errorSC = 2;
      goto cleanUp_main_sec04;

   cleanUp_main_sec04:;
      if(kmerCountsHeapAryST)
         freeHeapAry_kmerCnt(
            kmerCountsHeapAryST,
            numberSequencesSI
         );
      kmerCountsHeapAryST = 0;

      return errorSC;
} /*main*/
```

# get number matching kmers in structures

After making you kmerCnt structure and adding the kmers
  in you can then convert sequences to kmer arrays with
  the `ntToKmerAry_kmerCnt()` function. This kmer array
  can then be compared to a kmerCnt structure using
  the `get_kmerCnt()` function.

The `ntTokmerAry_kmerCnt()` function converts a sequence
  stored in a seqST structure to  an array of kmers.

- Input:
  1. seqST structure pointer with sequence to convert
  2. length of one kmer (should be same length used in
     making the  kmerCnt structure (input 2
     of `setup_kmerCnt` and `faToKmerCnt_kmerCnt`)
  3. signed int array to hold the kmer ids
     - size needs to be (kmer length * kmer length)
  4. signed int array to hold number of times a kmer was
     repeated
     - size needs to be (kmer length * kmer length)
- Output:
  - Returns number of kmers in the sequence
  - modifies input 3 and 4 to have kmer id's and counts

The `get_kmerCnt()` function compares the kmer counts
  from `ntToKmerAry_kmerCnt()` (input 3 and 4) to a
  kmerCnt structure. The return value is the number of
  kmers. The sign (positive or negative) of the returned
  count tells if the best match was forward or reverse
  complement.

- Input:
  1. kmerCnt structure pointer with reference to compare
  2. kmerArySI with kmer ids of the sequence
    - input 3 from `ntToKmerAry_kmerCnt()`
  3. cntArySI with the number of times each kmer was
     duplicated
    - input 4 from `ntToKmerAry_kmerCnt()`
- Output:
  - Returns positive number of matching kmers if the
    forward sequence had the most matches
  - Returns negative number of matching kmers if the
    reverse complement sequence had the most matches

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "seqST.h"
#include "kmerCnt.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec01: variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int *kmerHeapArySI = 0;
   signed int *countHeapArySI = 0;
   signed int kmerLenSI = 7; /*using 7mers*/

   signed int seqOnSI = 0;
   signed int numMatchesSI = 0;

   signed int errorSI = 0;
   struct seqST sequenceStackST;
   struct kmerCnt kmerCountsStackST;
   FILE *sequenceFILE = fopen("sequences.fasta", "r");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec02: initialize, setup, and get sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_kmerCnt(&kmerCountsStackST);
   init_seqST(&sequenceStackST);

   if(setup_kmerCnt(&kmerCountsStackST, kmerLenSI))
      goto memoryError_main_sec04;

   errorSI = getFa_seqST(sequenceFILE, &sequenceStackST);

   if(! errorSI)
      ; /*no errors*/
   else if(errorSI == def_memErr_seqST)
      goto memoryError_main_sec04;
   else
      goto fileError_main_sec04;

   /*initialize kmer arrays*/
   kmerHeapArySI =
      calloc(kmerLenSI * kmerLenSI, sizeof(signed int));
   if(! kmerHeapArySI)
      goto memoryError_main_sec04;

   countHeapArySI =
      calloc(kmerLenSI * kmerLenSI, sizeof(signed int));
   if(! countHeapArySI)
      goto memoryError_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec03: compare kmers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSI = getFa_seqST(sequenceFILE, &sequenceStackST);
   seqOnSI = 1;

   while(! errorSI)
   { /*Loop: compare rest of sequences to first sequence*/
     ++seqOnSI;

     ntToKmerAry_kmerCnt(
        &sequenceStackST, /*has sequence to convert*/
        kmerLenSI,        /*length of one kmer*/
        kmerHeapArySI,
        countHeapArySI
     ); /*convert seqeunces*/

     numMatchesSI =
        get_kmerCnt(
           &kmerCountsStackST,
           kmerHeapArySI,
           countHeapArySI
        );

     printf(
      "sequence %i has %i matching kmers to sequence 1\n",
      seqOnSI,
      numMatchesSI
     );

     errorSI = getFa_seqST(sequenceFILE,&sequenceStackST);
   } /*Loop: compare rest of sequences to first sequence*/

   if(errorSI == def_EOF_seqST)
      ; /*no errors*/
   else if(errorSI == def_memErr_seqST)
      goto memoryError_main_sec04;
   else
      goto fileError_main_sec04;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec04: return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSI = 0; /*no error*/
   goto cleanUp_main_sec04;

   memoryError_main_sec04:;
      fprintf(stderr, "ran out of memory\n");
      errorSI = 1;
      goto cleanUp_main_sec04;

   fileError_main_sec04:;
      fprintf(stderr, "error reading sequence.fasta\n");
      errorSI = 2;
      goto cleanUp_main_sec04;

   cleanUp_main_sec04:;
      freeStack_kmerCnt(&kmerCountsStackST);
      freeStack_seqST(&sequenceStackST);

      if(kmerHeapArySI)
         free(kmerHeapArySI);
      kmerHeapArySI = 0;

      if(countHeapArySI)
         free(countHeapArySI);
      countHeapArySI = 0;

      return errorSI;
} /*main*/

# clean up

When finished, remember to free you kmerCnt structure
  with the `freeStack_kmerCnt()` for stack allocated
  structures and `freeHeap_kmerCnt()` for heap allocated
  strucutres. Both functions take the pointer to the
  kmerCnt structure to free. Just remember to set your
  pointer to 0/null after using `freeHeap_kmerCnt()`.

For kmerCnt arrays you can use
  the `freeStackAry_kmerCnt()` to free stack allocated
  kmerCnt arrays or `freeHeapAry_kmerCnt()` to free heap
  allocated kmerCnt arrays. The input is a pionter to
  the kmerCnt array to free and the number (index 1) of
  structures in the array. For `freeHeapAry_kmerCnt()`
  remember to set your pointer to 0/null after calling.

See the examples in the previous sections for use of the
  free functions.
