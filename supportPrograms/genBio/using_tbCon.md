# Goal

Describe how to use the tbCon library to build a reference
  based consensus from a sam.

Files: tbCon.h tbCon.c tbConDefs.h

- Depenedencies:
  - samEntry
  - base10Str from genLib
  - ulCp from genLib
  - numToStr from genLib
  - tbConDefs.h
  - endLine.h from genLib
  - genMath.h from genLib

- Hidden dependencies:
  - fileFun from genLib
  - ntTo5Bit.h

The tbCon library using the samEntry structure, so you
  should be familar with the samEntry library (see
  using\_samEntry.md).

The tbConDefs.h file has the default settings in
  predefied error values.

# Overview

The tbCon system converts reads in a sam file into a
  consensus. However, indstead of aligning insertions, it
  requires an exact match. This results in poor handeling
  of insertions in ONT reads.

There are two structures to be familar with for praticle
  use. Then one structure if your really trying to
  understand how it all works.

- Two critical structures:
  1. set\_tbCon holds the settings for the consensus
  2. conNt\_tbCon holds a single base in the consensus
     - most of the time you will be passing this to
       functions, not working with it directly

# The flow

1. create a set\_tbCon structure and initialize to
   default settings
2. add reads in sam file to a conNt\_tbCon structure
   - functions will automatically do this for you, so
     you just create a conNt\_tbCon structure pointer and
     pass in the pointer and samEnty structure with the
     read to the funcction
3. collapse the consensus to get a samEntry structure with
   the consensus
   - and/or you can print the variants to a file
4. when finished with the consensus, free all structures

# Structures

- Structures:
  1. set\_tbCon holds the settings for the consensus
  2. conNt\_tbCon holds a single base in the consensus
     - most of the time you will be passing this to
       functions, not working with it directly
  3. ins\_tbCon is a single insertion event
     - this is included in conNt\_tbCon

## 01 settings

The set\_tbCon structure holds the settings for the
  consensus.

- Variables for with settings adding reads to the
  consensus
  1. minQSI: holds the  minimum q-score to add a single
     nucleotide to the consensus
     - does not apply to insertions
  2. minInsQSI: holds the minimum q-score to add a single
     nucleotide from an insertion to the consensus
     - minQSI for insertions
  3. minMapqUC: minimum mapping quality to add a read to
     the consensus
- Variables for with settings for collapsing the consensus
  1. maskSC: character to mask low read depth or low
     percentage support positions with
  2. minDepthSI: minimum read depth to keep a nucleotide
     in the collapsed consensus
     - for non-fragmentation; this masks the position
     - for a fragmented consensus (samEntry array) this
       breaks the consensus into fragments at low depths
  3. minPercSnpF: (0 to 1) minimum support needed to not
     mask a nucleotide in the collapsed consensus
  4. minPercDelF: (0 to 1) minimum support needed to not
     mask (remove) a deletion in the collapsed consensus
  5. minPercInsF: (0 to 1) minimum support needed to not
     an insertion in the collapsed consensus
  6. minLenSI: this is for the framgment consensus
     collapse, were an array of samEntry structures are
     returned
     - minimum length needed to keep a collapsed consensus
- Variables for printing variations to a file
  1. minPrintDepthSI: minimum variant depth to print
     the variant to a file
  2. printMinSupSnpF: minimum percent support to print a
     nucleotide variant (match/snp) to a file
  3. printMinSupInsF: minimum percent support to print a
     insertion variant to a file
  4. printMinSupDelF: minimum percent support to print a
     deletion variant to a file

### settings functions

There are only two functions for set\_tbCon. The first
  is `init_set_tbCon()`, which initializes all values to
  the defaults (only ever call once). The input is a
  pointer to the set\_tbCon structure to initialize.

The second is `freeStack_set_tbCon()`, which frees all
  variables in the set\_tbCon structure. The input is a
  pointer to the set\_tbCon structure to initialize.
  Currently, there is no memory stored, so all this does
  is call `init_set_tbCon()`. However, it is good pratice
  to use this, since unlikley future cases could could
  have memory to free.

```
#include "tbCon.h";

int
main(
){
   struct set_tbCon settingsStackST = 0;
   init_set_tbCon(&settingsStackST);

   /*do stuff here*/
   freeStack_set_tbCon(insHeapST);
   return 0;
}
```

## 02 nucleotide structure

The conNt\_tbCon holds a single nucleotide in a consensus.
  Again, you often will not have to do much directly with
  this other than pass it around to functions.

- Variables:
  1. numASI: number of reads that supported an A for this
     position
  2. numCSI: number of reads that supported a C for this
     position
  3. numGSI: number of reads that supported a G for this
     position
  4. numTSI: number of reads that supported a T for this
     position
  5. numDelSI: number of reads that supported a deletion
     for this position
  6. numInsSI: number of insertion events (1 or more
     nucleotides) at this position
     - tbCon always puts insertion on the nucleotide
       before the insertion
  7. insList: linked list of ins\_tbCon structures with
     the insertions (see next structure)
  8. totlaNtSI: number of reads that had a nucleotide or
     deletion at this position
  9. ntKeptSI: number of reads with nucleotides that
     passed the quality controll thresholds
     - also includes the deletion counts, which have no
       method to check quality

### work with the nucleotide structure

You can make a conNt\_tbCon structure using
  the `mk_conNt_tbCon()`. The input is nothing and the
  return value is a heap allocated conNt\_tbCon strucutre
  pointer. It will return 0 for memory errors.

You can then free your conNt\_tbCon structure using
  the `free_conNt_tbCon()` function. The input is the
  conNt\_tbCon structure to free. After freeing, remember
  to set you conNt\_tbCon structure to 0/null.

```
#include "tbCon.h";

int
main(
){
   struct conNt_tbCon *ntHeapST = 0;
   ntHeapST = mk_conNt_tbCon();
   if(! ntHeapST)
      /*deal with memory error*/
   
   /*do stuff here*/
   freeHeap_conNt_tbCon(insHeapST);
   ntHeapST = 0;

   return 0;
}
```

You can reset the values in a conNt\_tbCon structure back
  to nothing with the `blank_conNt_tbCon()`. This will set
  the variables to 0, but will do nothing with the
  ins\_tbCon structure.

`blank_conNt_tbCon(ntHeapST);`

To free the ins\_tbCon structures and set all values to
  0 in a conNt\_tbCon structure use
  the `freeStack_conNt_tbCon()` function. This function is
  also used to free stack allocated conNt\_tbCon
  structures.

`freeStack_conNt_tbcon(ntHeapST)`;

If you need to initialize a conNt\_tbCon structure, such
  as when you make an array, you can use
  the `init_conNt_tbCon()` function.

```
#include "tbCon.h";

int
main(
){
   signed int siNt = 0;
   struct conNt_tbCon ntArrayStackST[100];

   for(siNt = 0; siNt < 100; ++siNt)
      init_conNt_tbCon(&ntArrayStackST[siNt]);

   /*do stuff here*/

   for(siNt = 0; siNt < 100; ++siNt)
      freeStack_conNt_tbCon(&ntArrayStackST[siNt]);

   return 0;
}
```

Finally if you made a heap array of conNt\_tbCon
  structures, you can use the `freeHeapAry_conNt_tbcOn()`
  to free the array. The input is a pointer to the first
  conNtAryST in the array to free and the number of
  structures in the array.

```
#include "tbCon.h";

int
main(
){
   signed int siNt = 0;
   struct conNt_tbCon *ntArrayHeapST = 0;

   ntArrayHeapST =
      malloc(100 * sizeof(struct conNt_tbCon));
   if(! ntArrayHeapST)
      /*deal with memory error*/

   for(siNt = 0; siNt < 100; ++siNt)
      init_conNt_tbCon(&ntArrayStackST[siNt]);

   /*do stuff here*/

   freeHeapAry_conNt_tbCon(ntArrayHeapST, 100);
   ntArrayHeapST = 0;
      /*free(ntArrayHeapST) was called by freeHeapAry*/

   return 0;
}
```

## 03 insertion structure

The ins\_tbCon structure holds a single insertion event
  (one or more inserted nucleotides). As a rule, you
  should never have to dig this deep unless you are just
  using my structures or doing additinal analysis.

- Variables:
  1. insStr: c-string with the insertion sequence
  2. lenInsSI: length of the insertion
  3. numHitsSI: number of reads exactly matching this
     insertion
  4. nextIns: next insertion (ins\_tbCon struct) in the
     insertion list

### make insertion structure

You can make an ins\_tbCon structure using
  the `mkIns_tbCon()` function. There are not inputs and
  it returns 0 for memory errors or a pointer to a heap
  allocated ins\_tbcon structure.

You can free your ins\_tbCon structure with
  the `freeHeap_ins_tbCon()`. Just remember to set your
  pointer to 0/null when after freeing.

```
#include "tbCon.h";

int
main(
){
   struct ins_tbCon *insHeapST = 0;
   insHeapST = mkIns_tbCon();
   if(! insHeapST)
      /*deal with memory error*/
   
   /*do stuff here*/
   freeHeap_ins_tbCon(insHeapST);
   insHeapST = 0;
}
```

# using tbCon

There are four non-structure specific functions that can
  be used for tbCon.

1. addRead\_tbCon: adds a read to a conNt\_tbCon structure
   array
   - this function will always be used, unless you are
     doing your own system
   - this function will make the conNt\_tbCon array for
     you
2. collapse\_tbCon: collapases (makes a consensus) a
   conNt\_tbCon array into an array of samEntry structures
   with the consensus fragments
   - breaks up consnesus at points with low read depth
3. noFragCollapse\_tbCon: collapases a conensus and fills
   in low read depth positions with you masking character
   - this is like collapse\_tbCon, except output is only
     one samEntry structure instead of an array
4. pvar\_tbCon: prints the possible variants for each
   position in a conNt\_tbCon array
   - this does not generate a consensus and is the
     `-out-tsv` option for `tbCon`

## 1 add reads

The first step in building a consensus or getting variants
  is to add the reads to your consensus. You can do this
  with the `addRead_tbCon()` function.

- Input:
  1. samEntry structure pointer with read to add
  2. conNt\_tbCon structure pointer to add read to
     - this will auto make or expand the input when needed
  3. unsigned int pointer with the current length of the
     conNt\_tbCon array
     - again this is auto filled
  4. set\_tbCon strucuture pointer with the settings for
     adding reads (minimum Q-scores and mappng quality)
- Output:
  1. Modifies: conNt\_tbCon structure array (input 2) to
     the correct size
  2. Modifies the unsigned int (input 3) to have the
     expanded size of the conNt\_tbCon structure array
     (input 2)
  3. Returns:
     - 0 for no errors
     - def\_header\_tbConDefs if samEntry structure had
       a sam file header (not a read)
     - def\_noMap\_tbConDefs if read was unmapped
     - def\_noSeq\_tbConDefs if read had no sequence
     - def\_lowMapq\_tbConDefs if read had to low of a
       mapping quality
     - def\_memErr\_tbConDefs for memory errors

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "tbCon.h"
#include "tbConDefs.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   struct conNt_tbCon *conHeapSTAry = 0;
   unsigned int numBasesUI = 0;
   struct set_tbCon settingsStackST;
   struct samEntry readStackST;
   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   init_set_tbCon(&settingsStackST);
   init_samEntry(&readStackST);
   if( setup_samEntry(&readStackST) )
      goto memErr_main_sec05;
   
   samFILE = fopen("reads.sam", "r");
   if(! samFILE)
      goto fileErr_main_sec05;
   
   errorSC = get_samEntry(&readStackST, samFILE);
   
   if(! errorSC)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec05;
   else
      goto fileErr_main_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: add reads to consensus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   while(! errorSC)
   { /*Loop: add reads to the consnesus*/
     
      /*for simipilicity I am assuming everything mapps
      `  to the same reference
      */
      errorSC =
         addRead_tbCon(
            &readStackST,
            &conHeapSTAry, /*function auto gets memory*/
            &numBasesUI,   /*function auto gets this*/
            &settingsStackST     
         );

      if(errorSC == def_memErr_tbConDefs)
         goto memErr_main_sec05;

      errorSC = get_samEntry(&readStackST, samFILE);
   } /*Loop: add reads to the consnesus*/

   if(errorSC != def_EOF_samEntry)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec05;
   else
      goto fileErr_main_sec05;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: finish consnesus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*will cover these functions next*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto ret_main_sec05;

   memErr_main_sec05:;
      errorSC = 1;
      fprintf(stderr, "memory error\n");
      goto ret_main_sec05;

   fileErr_main_sec05:;
      errorSC = 2;
      fprintf(stderr, "file error with reads.sam\n");
      goto ret_main_sec05;

   ret_main_sec05:;
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
      
      freeStack_set_tbCon(&settingsStackST);
      freeStack_samEntry(&readStackST);
      freeHeapAry_conNt_tbCon(conHeapSTAry, numBasesUI);

      return errorSC;
```

## 2 collapse consensus into fragments

For collapsing the consensus you can collapse it into
  fragments, were each position with low read depth is
  broken into a separate consensus. The function to do
  this is `collapse_tbCon()`.

- Input
  1. an array of conNt\_tbCon structures to collapse
  2. number if items in the conNt\_tbCon array (input 1)
  3. signed int pointer to get the number of fragments
     that were made
  4. c-string with the name of the reference used
  5. set\_tbCon structure with settings for the collapse
     - minLenSI, minDepthSI, minPercSnpF, minPercDelF,
       minPercInsF, maskSC
  6. signed character to hold the error type
- Output:
  1. Modifies: singed int pionter (input 3) to have the
     number of consensuses made
  2. Modifies: singed character pionter (input 6) to have
     the returned error type
     - 0 for no errors
     - def\_noSeq\_tbConDefs for could not collapse the
       consensues
     - def\_memErr\_tbConDefs for memory errors
  3. Returns:
     - array of samEntry structures with the consensus
     - 0 if had an error

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "tbCon.h"
#include "tbConDefs.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   struct conNt_tbCon *conHeapSTAry = 0;
   unsigned int numBasesUI = 0;
   struct set_tbCon settingsStackST;
   struct samEntry readStackST;
   FILE *samFILE = 0;

   struct samEntry *consensusesHeapST = 0;
   signed int numConsSI = 0;
   signed int siCon = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   init_set_tbCon(&settingsStackST);
   init_samEntry(&readStackST);
   if( setup_samEntry(&readStackST) )
      goto memErr_main_sec05;
   
   samFILE = fopen("reads.sam", "r");
   if(! samFILE)
      goto fileErr_main_sec05;
   
   errorSC = get_samEntry(&readStackST, samFILE);
   
   if(! errorSC)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec05;
   else
      goto fileErr_main_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: add reads to consensus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   while(! errorSC)
   { /*Loop: add reads to the consnesus*/
     
      /*for simipilicity I am assuming everything mapps
      `  to the same reference
      */
      errorSC =
         addRead_tbCon(
            &readStackST,
            &conHeapSTAry, /*function auto gets memory*/
            &numBasesUI,   /*function auto gets this*/
            &settingsStackST     
         );

      if(errorSC == def_memErr_tbConDefs)
         goto memErr_main_sec05;

      errorSC = get_samEntry(&readStackST, samFILE);
   } /*Loop: add reads to the consnesus*/

   if(errorSC != def_EOF_samEntry)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec05;
   else
      goto fileErr_main_sec05;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: finish consnesus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   consensusesHeapST =
      collapse_tbCon
         conHeapSTAry, /*to collapse; from addRead_tbCon*/
         numBasesUI,   /*length of conHeapSTAry*/
         &numConsSI,
         (signed char *) "reference_name",
         &settingsStackST, /*settings*/
         &errorSC      /*gets errors*/
      );

   if(! errorSC)
      ; /*no errors*/
   else if(errorSC == def_memErr_tbConDefs)
      goto memErr_main_sec05;
   else
      goto noConErr_main_sec05;

   for(siCon = 0; siCon < numConsSI; ++siCon)
   { /*Loop: print out consensus*/
      if(consensusesHeapST[siCon].qryIdStr[0] != '\0')
         p_samEntry(&consensusHeapST[siCon], 0, stdout);
         /*this is needed to avoid blank entries; the
         `  fragmentation system is not perfect
         */
   } /*Loop: print out consensus*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto ret_main_sec05;

   memErr_main_sec05:;
      errorSC = 1;
      fprintf(stderr, "memory error\n");
      goto ret_main_sec05;

   fileErr_main_sec05:;
      errorSC = 2;
      fprintf(stderr, "file error with reads.sam\n");
      goto ret_main_sec05;

   noConErr_main_sec05:;
      errorSC = 3;
      fprintf(stderr, "could not build consensus\n");
      goto ret_main_sec05;

   ret_main_sec05:;
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
      
      if(consensusesHeapST)
      { /*If: have to free the consensus fragments*/
         for(siCon = 0; siCon < numConsSI; ++siCon)
            freeStack_samEntry(&consensusesHeapST[siCon]);
         free(consensusesHeapST);
      } /*If: have to free the consensus fragments*/
      consensusesHeapST = 0;

      freeStack_set_tbCon(&settingsStackST);
      freeStack_samEntry(&readStackST);
      freeHeapAry_conNt_tbCon(conHeapSTAry, numBasesUI);

      return errorSC;
```

## 3 collapse consensus

If you want a single consensus sequence, you can collapse
  it using the `noFragCollapse_tbCon()` function.

- Input
  1. an array of conNt\_tbCon structures to collapse
  2. number if items in the conNt\_tbCon array (input 1)
  4. c-string with the name of the reference used
  3. unsigned int pointer to get the number of bases that
     that were masked
  5. set\_tbCon structure with settings for the collapse
     - minDepthSI, minPercSnpF, minPercDelF, minPercInsF,
       maskSC
  6. signed character to hold the error type
- Output:
  1. Modifies: singed int pionter (input 3) to have the
     number of consensuses made
  2. Modifies: singed character pionter (input 6) to have
     the returned error type
     - 0 for no errors
     - def\_noSeq\_tbConDefs for could not collapse the
       consensues
     - def\_memErr\_tbConDefs for memory errors
  3. Returns:
     - a heap allocated samEntry structure with the
       consensus
     - 0 if had an error

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "tbCon.h"
#include "tbConDefs.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   struct conNt_tbCon *conHeapSTAry = 0;
   unsigned int numBasesUI = 0;
   struct set_tbCon settingsStackST;
   struct samEntry readStackST;
   FILE *samFILE = 0;

   struct samEntry *conHeapST = 0;
   unsigned int numMaskUI = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   init_set_tbCon(&settingsStackST);
   init_samEntry(&readStackST);
   if( setup_samEntry(&readStackST) )
      goto memErr_main_sec05;
   
   samFILE = fopen("reads.sam", "r");
   if(! samFILE)
      goto fileErr_main_sec05;
   
   errorSC = get_samEntry(&readStackST, samFILE);
   
   if(! errorSC)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec05;
   else
      goto fileErr_main_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: add reads to consensus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   while(! errorSC)
   { /*Loop: add reads to the consnesus*/
     
      /*for simipilicity I am assuming everything mapps
      `  to the same reference
      */
      errorSC =
         addRead_tbCon(
            &readStackST,
            &conHeapSTAry, /*function auto gets memory*/
            &numBasesUI,   /*function auto gets this*/
            &settingsStackST     
         );

      if(errorSC == def_memErr_tbConDefs)
         goto memErr_main_sec05;

      errorSC = get_samEntry(&readStackST, samFILE);
   } /*Loop: add reads to the consnesus*/

   if(errorSC != def_EOF_samEntry)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec05;
   else
      goto fileErr_main_sec05;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: finish consnesus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   conHeapST =
      noFragCollapse_tbCon
         conHeapSTAry, /*to collapse; from addRead_tbCon*/
         numBasesUI,   /*length of conHeapSTAry*/
         (signed char *) "reference_name",
         &numMaskUI,
         &settingsStackST, /*settings*/
         &errorSC      /*gets errors*/
      );

   if(! errorSC)
      ; /*no errors*/
   else if(errorSC == def_memErr_tbConDefs)
      goto memErr_main_sec05;
   else
      goto noConErr_main_sec05;

   p_samEntry(conHeapST, 0, stdout);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto ret_main_sec05;

   memErr_main_sec05:;
      errorSC = 1;
      fprintf(stderr, "memory error\n");
      goto ret_main_sec05;

   fileErr_main_sec05:;
      errorSC = 2;
      fprintf(stderr, "file error with reads.sam\n");
      goto ret_main_sec05;

   noConErr_main_sec05:;
      errorSC = 3;
      fprintf(stderr, "could not build consensus\n");
      goto ret_main_sec05;

   ret_main_sec05:;
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
      
      if(conHeapST)
         freeHeap_samEntry(conHeapST);
      conHeapST = 0;

      freeStack_set_tbCon(&settingsStackST);
      freeStack_samEntry(&readStackST);
      freeHeapAry_conNt_tbCon(conHeapSTAry, numBasesUI);

      return errorSC;
```

## printing variants in a consensus

If you want a list of possible variants for the conensus
  you can use the `pvar_tbCon()` function. This will
  report the variant, total read depth, variant read
  depth, and the percent support for the varaint.

- Input
  1. an array of conNt\_tbCon structures to get variants
     from
  2. number if items in the conNt\_tbCon array (input 1)
  4. c-string with the name of the reference used
  5. set\_tbCon structure with settings for the collapse
     - minPrintDepthSI, minPrintSupSnpF, minPrintSupInsF,
       and minPrintSupDelF
  6. c-string with file name to print to or 0 for stdout
- Output:
  - Returns:
     - 0 for no problems
     - def\_fileErr\_tbConDefs for file errors
     - def\_noMap\_tbConDefs if no conNt\_tbCon array
       input (or is 0 (null))

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "tbCon.h"
#include "tbConDefs.h"

int
main(
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01: variable declrations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errorSC = 0;
   struct conNt_tbCon *conHeapSTAry = 0;
   unsigned int numBasesUI = 0;
   struct set_tbCon settingsStackST;
   struct samEntry readStackST;
   FILE *samFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02: initialize and setup structures
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   init_set_tbCon(&settingsStackST);
   init_samEntry(&readStackST);
   if( setup_samEntry(&readStackST) )
      goto memErr_main_sec05;
   
   samFILE = fopen("reads.sam", "r");
   if(! samFILE)
      goto fileErr_main_sec05;
   
   errorSC = get_samEntry(&readStackST, samFILE);
   
   if(! errorSC)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec05;
   else
      goto fileErr_main_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03: add reads to consensus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   while(! errorSC)
   { /*Loop: add reads to the consnesus*/
     
      /*for simipilicity I am assuming everything mapps
      `  to the same reference
      */
      errorSC =
         addRead_tbCon(
            &readStackST,
            &conHeapSTAry, /*function auto gets memory*/
            &numBasesUI,   /*function auto gets this*/
            &settingsStackST     
         );

      if(errorSC == def_memErr_tbConDefs)
         goto memErr_main_sec05;

      errorSC = get_samEntry(&readStackST, samFILE);
   } /*Loop: add reads to the consnesus*/

   if(errorSC != def_EOF_samEntry)
      ;
   else if(errorSC == def_memErr_samEntry)
      goto memErr_main_sec05;
   else
      goto fileErr_main_sec05;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04: finish consnesus
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC =
      pvar_tbCon
         conHeapSTAry, /*to collapse; from addRead_tbCon*/
         numBasesUI,   /*length of conHeapSTAry*/
         (signed char *) "reference_name",
         &settingsStackST, /*settings*/
         0             /*outputing to stdout*/
      );

   if(! errorSC)
      ; /*no errors*/
   else if(errorSC == def_fileErr_tbConDefs)
      goto fileErr_main_sec05;
   else
      goto noConErr_main_sec05;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec05: clean up and return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errorSC = 0;
   goto ret_main_sec05;

   memErr_main_sec05:;
      errorSC = 1;
      fprintf(stderr, "memory error\n");
      goto ret_main_sec05;

   fileErr_main_sec05:;
      errorSC = 2;
      fprintf(stderr, "file error with reads.sam\n");
      goto ret_main_sec05;

   noConErr_main_sec05:;
      errorSC = 3;
      fprintf(stderr, "could not build consensus\n");
      goto ret_main_sec05;

   ret_main_sec05:;
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
      
      freeStack_set_tbCon(&settingsStackST);
      freeStack_samEntry(&readStackST);
      freeHeapAry_conNt_tbCon(conHeapSTAry, numBasesUI);

      return errorSC;
```
