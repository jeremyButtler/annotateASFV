# Goal:

Describe how to use the alnSet library to set up for the
  Neeldeman, Waterman, or Hirschberg aligners in bioTools.

Files: alnSet.h and alnSet.c

- Dependencies:
  - alnDefs.h
  - base10str from genLib
  - ulCp from genLib
  - fileFun from genLib
  - endLine.h from genLib

# alnSet overview:

The alnSet library has the structure and supporting
  functions to store the alignment matrix, scoring
  penalties. The functions will allow you to change
  scores, if something is a match, and convert sequences
  to the correct index for alignment.

# Structure:

The alnSet structure stores the alignment scoring matrix
  and other alignment variables.

- Variables:
  - minScoreSL: minimum score to keep alignment (Waterman)
  - gapSS: gap opening penalty (should be negative)
  - extendSS: gap extension penalty (should be negative)
  - delArySS: used to get scores for extending deletions,
              or adding deletion after an insertion, snp
              or match
  - insArySS: used to get scores for extending insertions,
              or adding insetion after a deletion, snp
              or match
  - socreMatrixSS: scoring matrix for alignment
  - matchMatrixSS: decides if two nucleotides are matches
                   or SNP

The delArySS and insArrySS use index's from the movement
  defined variables in alnDef. This makes since, since
  in an aligment you are keeping track of the last
  direction (mutation type).

- These movement (mutation type) variables are
  - def\_mvStop\_alnDef: stop movement (for waterman)
  - def\_mvDel\_alnDef: deletion movement
  - def\_mvIns\_alnDef: insertion movement
  - def\_mvSnp\_alnDef: mismatch/match movement

```
unsigned char firstDirUC = def_mvDel_alnDefs;
unsigned char secDirUC = def_mvIns_alnDefs;
struct alnSet alnStackST;

init_alnSet(&alnStackST);

printf("score for deletion when first movement");
printf(
   " is a deletion and second is insertion %i\n",
   alnStackST.delArySS[firstDirUC][secDirUC]
);

printf("score for insertion when first movement");
printf(
   " is a deletion and second is insertion %i\n",
   alnStackST.insArySS[firstDirUC][secDirUC]
);
```

# Functions:

## initialize and freeing

To use an alnSet structure you must first initialize it
  to the defaults in the alnDefs.h file. You can do this
  using the `init_alnSet()` function.

- Input:
  1. alnSet structure pointer to initialize

When finished you can free the alnSet structure using
  the `freeStack_alnSet()` function to free the variables
  in an alnSet structure (does nothing, here for future
  use).

You can also use the `freeHeap_alnSet()` function to free
  a heap allocated alnSet structure. After calling you
  should set your pointer to null/0.

```
struct alnSet alnStackST;
init_alnSet(&alnStackST);

/*code goes here*/

freeStack_alnSet(&alnStackST);
```

## set matrix scores and matches

For setting scores and matches, you can set the scores
  from a file or you can set each score/match one by one.

### setting scores individually

You can set each base pair score individually in the
  alnSet structures scoring matrix using the
  the `setScore_alnSet()` macro.

- Input:
  1. first (query) nucleotide in pair
     - can be as character or alnSet code
  2. second (reference) nucleotide in pair
     - can be as character or alnSet code
  3. new score to assign to the base pair
  4. alnSet structure pionter to add the score to

**Example**

```
struct alnSet alnStackST;
init_alnSet(&alnStackST);
setScore_alnSet('A', 'G', -2, &alnStackST);
```

### setting matches individually

You can set each base pair match individually in the
  alnSet structures match matrix using the
  the `setScore_alnSet()` macro.

- Input:
  1. first (query) nucleotide in pair
     - can be as character or alnSet code
  2. second (reference) nucleotide in pair
     - can be as character or alnSet code
  3. match to assign
     - 0 for not match
     - 1 for match
     - def\_anonymous\_alnDefs for anonymous bases
     - def\_anonMatch\_alnDefs for anonymous bases that
       is a match
  4. alnSet structure pionter to add the score to

**Example**

```
struct alnSet alnStackST;
init_alnSet(&alnStackST);
setMatch_alnSet('A', 'A', 1, &alnStackST); /*match*/
setMatch_alnSet('A', 'G', 0, &alnStackST); /*mismatch*/

setMatch_alnSet('A','W',def_anonMatch_alnDefs,&alnStackST);
   /*match for an anonymous base*/
setMatch_alnSet('A','K',def_anonymous_alnDefs,&alnStackST);
   /*mismatch for an anonymous base*/
```

### getting scores from a file

You can read in base pair scores for the alignment with
  the `readScoreFile_alnSet()` function.

- Input:
  1. alnSet structure pointer to load the scoring matrix
     into
  2. FILE pointer with file to get scores from
- Output:
  - Returns: 0 for no errors
  - Returns: line number in file had error on

The file format for the scoring file has one line per
  comparision. Each line in the file must be blank, be a
  comment (start with a `//`), or
  have the score for the pair `<query> <reference> score`.

For example, to set tranvsions to -4 and transitions to
  -1 for `A` I could do.

```
// This is  a comment

A A 5   // match (this is another comment)
A T -4  // A to T is a transversion
A U -4  // A to T is a transversion
A C -4  // A to C is a transversion
A G -1  // A to G is a transition
```

**Example**

```
#include <stdio.h>
#include "alnSet.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec01: variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned long errUL = 0;
   struct alnSet alnStackST;
   FILE *scoreFILE = fopen("scoring-matrix.txt", "r");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec02: initialize and get scores from file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_alnSet(&alnStackST);
   
   if(! scoreFILE)
      goto fileErr_main_sec03;
   
   errUL = readScoreFile_alnSet(&alnStackST, scoreFILE);

   if(errUL)
   { /*If: had file error*/
      fprintf(stderr, "error on line %ul\n", errUL);
      goto fileErr_main_sec03;
   } /*If: had file error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec03: return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errUL = 0;
   goto ret_main_sec03;

   fileErr_main_sec03:;
      errUL = 1;
      goto ret_main_sec03;
   
   ret_main_sec03;
      freeStack_alnSet(&alnStackST);
   
      if(! scoreFILE) ;
      else if(scoreFILE == stdin) ;
      else if(scoreFILE == stdout) ;
      else if(scoreFILE == stderr) ;
      else
         fclose(scoreFILE);
      scoreFILE = 0;
   
      return errUL;
} /*main*/
```

### getting matches from a file

You can read in which base pairs are matches for the
  alignment with the `readMatchFile_alnSet()` function.

- Input:
  1. alnSet structure pointer to load the scoring matrix
     into
  2. FILE pointer with file to get matches from
- Output:
  - Returns: 0 for no errors
  - Returns: line number in file had error on


The file format for the matching file has one line per
  comparision. Each line in the file must be blank, be a
  comment (start with a `//`), or have if the pair
  is a match (1), mismatch (0), or anonymous
  (N) `<query> <reference> <code>`.

- Codes:
  - 1 for match
  - 0 for mismatch
  - N for anonymous base

For example, for the base A I might do.

```
// This is  a comment

A A 1   // match (this is another comment)
A T 0   // A to T is a mismatch
A U 0   // A to T is a mismatch
A C 0   // A to C is a mismatch
A G 0   // A to G is a mismatch
A S N   // S is an anonymous base
A N N   // N is an anonymous base
```

**Example**

```
#include <stdio.h>
#include "alnSet.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec01: variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned long errUL = 0;
   struct alnSet alnStackST;
   FILE *scoreFILE = fopen("scoring-matrix.txt", "r");

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec02: initialize and get scores from file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_alnSet(&alnStackST);
   
   if(! scoreFILE)
      goto fileErr_main_sec03;
   
   errUL = readMatchFile_alnSet(&alnStackST, scoreFILE);

   if(errUL)
   { /*If: had file error*/
      fprintf(stderr, "error on line %ul\n", errUL);
      goto fileErr_main_sec03;
   } /*If: had file error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Sec03: return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errUL = 0;
   goto ret_main_sec03;

   fileErr_main_sec03:;
      errUL = 1;
      goto ret_main_sec03;
   
   ret_main_sec03;
      freeStack_alnSet(&alnStackST);
   
      if(! scoreFILE) ;
      else if(scoreFILE == stdin) ;
      else if(scoreFILE == stdout) ;
      else if(scoreFILE == stderr) ;
      else
         fclose(scoreFILE);
      scoreFILE = 0;
   
      return errUL;
} /*main*/
```

## get matrix scores and matches

You can get the scores or match/mismatch call for bases
  using macros in the alnSet library.

### getting scores from the matrix

You can get the score for a base pair score from the
  alnSet structures scoring matrix using the
  the `getScore_alnSet()` macro.

- Input:
  1. first (query) nucleotide in pair
     - must be in lookup index
  2. second (reference) nucleotide in pair
     - must be in lookup index
  3. alnSet structure pionter to add the score to

**Example**

```
struct alnSet alnStackST;
init_alnSet(&alnStackST);

printf(
   "A and G have a score of: %i\n",
   getScore_alnSet(
      'A' & def_ntToCode_alnSet,
      'G' & def_ntToCode_alnSet,
       &alnStackST
   ) /*need ntToCode to convert A and G to the index for
     `  alnSet. There is a function to do this for
     `  c-strings
     */
);

freeStack_alnSet(&alnStackST);
```


### getting matches from the matrix

You can get if two bases are a match with the match matrix
  in the alnSet structu using the the `getMatch_alnSet()`
  macro.

- Input:
  1. first (query) nucleotide in pair
     - must be in lookup index
  2. second (reference) nucleotide in pair
     - must be in lookup index
  3. alnSet structure pionter to add the score to

**Example**

```
signed int matchSI = 0;
struct alnSet alnStackST;
init_alnSet(&alnStackST);

matchSI =
   getMatch_alnSet(
      'A' & def_ntToCode_alnSet,
      'G' & def_ntToCode_alnSet,
       &alnStackST
   ); /*need ntToCode to convert A and G to the index for
      `  alnSet. There is a function to do this for
      `  c-strings
      */

if(matchSI == 1)
   printf("A and G are a match\n"); /*never fires*/
else if(matchSI == 0)
   printf("A and G are a mismatch\n"); /*fires*/
else
   printf("A or G is an anonymous base\n");/*never fires*/

freeStack_alnSet(&alnStackST);
```

## set gap scores

You can change the gap opening and gap extension scores
  with the `changeGap_alnSet()` function. This will
  updated the gap open and gap extension scores in
  delArySS and insArySS. This will also set gapSS and
  extendSS.

- Input:
  1. alnSet structure pointer to add new score to
  2. new gap opening penalty (should be negative)
  3. new gap extension penalty (should be negative)

```
signed int matchSI = 0;
struct alnSet alnStackST;
init_alnSet(&alnStackST);

changeGap_alnSet(&alnStackST, gapSS, extendSS);

freeStack_alnSet(&alnStackST);
```

## convert sequences

Inorder to use the matrix's in an alnSet structure your
  need to convert the sequences to index. The alnSet
  library has a variable to convert nucleotides to the
  index format, a function to convert a sequence to
  index's, and a function to convert sequences back to
  nucleotides.

### convert nucleotide to index or index to nucleotide

You can convert a nucleotide to an index uinsing the
  def\_ntToCode\_alnSet variable to remove the higher
  bits. To convert back, you would add in 64
  (removes case).

**Example**

```
/*convert nucleotide to index*/
signed char indexSC = 'A' & def_ntToCode_alnSet;

/*to convert back to upper case*/
signed char nucleotides = indexSC | 64;
```

### convert nucleotid sequence to index sequence

You can convert an entire sequence to index's using the
  `seqToIndex_alnSet()` function. This modifies the input
  c-string to be in index format. This format is needed
  to use the score and match matrix in alnSet.

- Input:
  1. c-string with sequence to convert

You can convert the index sequence back to a nucleotide
  sequence using the `indexToSeq_alnSet()` function. This
  modifies the input c-string to be a uppercase nucleotide
  sequence.

- Input:
  1. c-string with sequence to convert

**Example**

```
signed char querySeqStr[128] = "ATGCCATGACCGA";
signed char referenceSeqStr[128] = "TTGCGATGTCGGA";
struct alnSet alnStackST;

init_alnSet(&alnStackST);

/*convert sequence to index's so can get socres*/
seqToIndex_alnSet(querySeqStr);
seqToIndex_alnSet(referenceSeqStr);

getScore_alnSet(
   querySeqStr[0],
   referenceSeqStr[0],
   &alnStackST
); /*get score for the first bases in both sequences*/

/*convert index's back to sequences*/
indexToSeq_alnSet(querySeqStr);
indexToSeq_alnSet(referenceSeqStr);

freeStack_alnSet(&alnStackST);
```

## misc

These functions do not fall in a catagory, but might be
  usefull.

### getting the maximum score for a sequence

You can get the maximum score possible for a sequence use
 the `maxScore_alnSet()` function.

- Input:
  1. c-string with sequence (as index) to get score for
     - use `seqToIndex_alnSet()` to convert to index
  2. first base to get score for (if unsure, use 0)
     - index 0
  3. length of sequence
     - index 1
  4. alnSet structure pointer with scoring matrix
- Output:
  - Returns: maximum score for sequence

**Example**

```
signed char seqStr[128] = "ATGCCATGACCGA";
struct alnSet alnStackST;

init_alnSet(&alnStackST);

/*convert sequence to index's so can get socres*/
seqToIndex_alnSet(seqStr);

printf(
   "max score for %s is %li",
   seqStr,
   maxScore_alnSet(seqStr, 0, 13, &alnStackST)
);
```

### reverse complementing sequences (as index)

You can get reverese complement a sequence in index format
 using the `revCmpIndex_alnSet()` function.

- Input:
  1. c-string with sequence (as index) to reverse
     complement
     - use `seqToIndex_alnSet()` to convert to index
  2. c-string to keep in sync with input 1
     - Example: quality score entry
     - use 0 for nothing
  3. length of sequence
     - index 1

**Example**

```
signed char seqStr[128] = "ATGCCATGACCGA";
print("sequence %s\n", seqStr);

/*you should never reverse complement this way (it is
`  inefficent), instead use seqST (see
`  ../genBio/using_seqST.md). This should only be used if
`  you plan to use the reverse complemented index's or
`  have already converted the sequences
*/
seqToIndex_alnSet(seqStr);
revCmpIndex_alnSet(seqStr, 0, 13);
indexToSeq_alnSet(seqStr);

print("reverse complement sequence %s\n", seqStr);
```
