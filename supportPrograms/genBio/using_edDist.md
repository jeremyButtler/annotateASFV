# Use:

Here to familarize your with the edDist library in genBio.

The edDist programs are used in finding edit distances
  between a read and a reference or between two reads
  mapped to the same reference (less accurate).

Files: edDist.c and edDist.h

- Dependencies:
  - ulCp from genLib
  - endLine.h from genLib
  - genMath.h from genLib
  - seqST from genBio
  - samEntry from genBio

- Hidden dependencies:
  - base10str from genLib
  - numToStr from genLib
  - fileFun from genLib
  - ntTo5Bit from genBio

You should be familar with the samEntry (see
  using_samEntry.md) and seqST (see using_seqST) before
  using edDist.

# edDist:

Get a modified edit distance that adjusts for indels
  being error prone (ignores small indels).

The results are often stored in a `res_edDist` structure.

General error types:

| Error Type              | Error Macro (name)   |
|:------------------------|:---------------------|
| file error              | def\_fileErr\_edDist |
| memory error (ran out)  | def\_memErr\_edDist  |
| no sequence in samEntry | def\_noSeq\_edDist   |
| read was not mapped     | def\_noMap\_edDist   |

## res\_edDist structure

This structure holds the stats for the edit distance
  function results. Some variables are only used for
  comparing two reads (less accurate).

- Result variables:
  - edDistSL: found edit distance
  - probUI: has the percent distance for the read
  - overlapUI: has the overlap shared between two reads
  - startUI: is the first base in the queyr that overlaps
             with the reference reads; comparing two reads
  - numIndelUI: number of indels counted
  - indelEventsUI: number of indel groups (events)
  - maxWinDisUI: maximum distance seen in all windows
  - minWinDisUI: minimum distance seen in all windows
  - avgWinDisUI: average distance seen between all windows
  - probMaxWinUI: maximum percent distance (like probUI)
                  seen in all the windows
             
## flow of a res\_edist structure

To setup a res\_edDist structure you must first make your
  structure `struct res\_edDist;`. Then you must
  initialize your `res\_edDist` structure. Finally, you
  can then free the `res\_edDist` structure at the end.

You can initialize with init\_res\_edDist (fun02). The
  input is a `res\_edDist` structure pointer. This
  function should only be called once.

You can free your `res_edDist` structure with
  `freeStack_res_edDist` (fun03) or `freeHeap_res_edDist`
  (fun04). FreeStack will free all variables inside the
  structure, while freeHeap will free all variables and
  the structure.
 
You can blank the structure with `blank_res_edDist`
  (fun01). However, both `readCmpDist_edDist`
  and `dist_edDist` will both call `blank_res_edDist`. So,
  only use the blank function for your own functions.

Here is an example for a stack workflow:

```
struct res_edDist distanceStackST;
init_res_edDist(&distanceStackST);

/*your code goes here*/

freeStack_res_edDist(&distanceStackST);
```

Here is an example for a heap workflow:

```
struct res_edDist *distanceHeapST = 0;
distanceHeapST = malloc(sizeof(struct res_edDist));
if(! distanceHeapST)
   /*handle memory error*/
init_res_edDist(distanceHeapST);

/*your code goes here*/

freeHeap_res_edDist(distanceHeapST);
```

## Finding edit distances

You can find the edit distance to a reference with
  `dist_edDist` or between two reads (less accurate) with
  `readCmpDist_edDist`. To get the percent distance use
  `percDist_edDist`.

### reference edit distance

The reference edit distance can be found with
  `dist_edDist` (fun07). It takes several inputs and
  returns the edit distance between the read (query) and
  reference sequence.

- Input:
  - samEntry structure pointer with the mapped read
  - seqST structure pointer with the reference sequence
    - this is to handle the `M` entries in the cigars,
      you can either use an eqx cigar, give a reference,
      or not count SNPs
      - you can get a eqx cigare with cigToEqx (see
        using\_cigToEqx.md)
      - in many cases your read mappers allow you to
        output eqx cigars (minimap2 is `--eqx`)
    - or 0 if 
  - minimum length to count an indel
    - any indel beneath this length is ingored
  - minimum quality score to count a base in the edit
    distance
  - size of one window
    - edDist counts distance by windows (ex 500 base
      windows) and total distance
  - `res_edDist` structure to store the extra results in
    - ex: window results, overlap, and indel counts

- Returns:
  - edit distance (0 or greater) for no errors
  - def\_noSeq\_edDist (is a negative number) if the
    read has no sequence
    - this would be a secondary alignment
  - def\_noMap\_edDist (is a negative number) if the
    read did not map to the reference genome

```
#ifdef PLAN9
   #include <u.h>
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "seqST.h"
#include "samEntry.h"
#include "edDist.h"

int
main(
){
   signed long distSL = 0;
   
   struct seqST refStackST
   struct res_edDist distanceStackST;
   struct samEntry readStackST;
   
   FILE *refFILE = fopen("reference.fa", "r");
   FILE *samFILE = fopen("reads.sam", "r");
   
   init_res_edDist(&distanceStackST);
   init_seqST(&refStackST);
   init_samEntry(&readStackST);
   
   if( setup_samEntry(&readStackST) )
      /*handel memory error*/
   
   if( getFa_seqST(refFILE, &refStackST) )
      /*handle the file error or memory error*/
   
   fclose(refFILE);
   refFILE = 0;
   
   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: read sam file*/
      distSL =
         dist_edDist(
            &readStackST,/*read to get distance for*/
            &refStackST, /*reference sequence mapped to*/
            5,         /*indels be at least 5 bases long*/
            7,         /*base has q-score of 7 or more*/
            500,       /*500 bases per window*/
            &distanceStackST
         );
   
      if(distSL == def_noSeq_edDist)
         continue; /*header or secondary alignment*/
      else if(distSL < 0)
         continue; /*unmapped read */
   
      printf(
         "distance for %s: %li\n",
         readStackST.qryIdStr,  /*read id*/
         distSL                 /*edit distance*/
      );
   } /*Loop: read sam file*/
   
   freeStack_res_edDist(&distanceStackST);
   freeStack_seqST(&refStackST);
   freeStack_samEntry(&readStackST);

   fclose(samFILE);
   samFILE = 0;

   return 0;
}
```

### getting read depths

You can get the depth of all reads to a single read in a
  sam file by using  `mkDepthProfile_edDist` (fun10). This
  function will reset the FILE location to the start of
  the sam file.

- Input:
  - samEntry structure pionter with the read to compare
    (or sequence) to compare all other reads to
  - minimum Q-score to count a single base in the depth
  - minimum overlap between read and reference (0 to 1) to
    count read in the read depth
  - `res_edDist` structure to store the read depths in
     - profile is stored in depthAryUI 
  - samEntry structure pointer used in reading each line
    of the file (avoid extra memory allocations)
  - FILE pointer to the sam file to get the read depths
    for

- Output:
  - 0 for no errors
  - def\_memErr\_edDist for no memory errors
  - def\_fileErr\_edDist for file errors 

Another similar function is `addReadToDepth_edDist`
  (fun09), were you add reads individually to the depth
  profile instead of scanning the entire file.

- Input:
  - samEntry structure pointer with reference sequence
    (ex other read) to compare to
  - samEntry pointer with read to add to the profile
  - minimum Q-score to count a single base in the depth
  - minimum overlap between read and reference (0 to 1) to
    count read in the read depth
  - `res_edDist` structure to store the read depths in
     - profile is stored in depthAryUI 
- Output:
  - 0 for success
  - 1 for failure

This example only shows you how to call the
  mkDepthProfile_edDist function. You would likley use
  this in finding read to read edit distances. See the
  next section for that example, but without the profile
  step.
  

```
#ifdef PLAN9
   #include <u.h>
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "edDist.h"

int
main(
){
   signed long distSL = 0;
   
   struct res_edDist distanceStackST;
   struct samEntry readStackST;
   struct samEntry refStackST
   
   FILE *samFILE = fopen("reads.sam", "r");
   
   init_res_edDist(&distanceStackST);
   init_samEntry(&refStackST);
   init_samEntry(&readStackST);
   
   if( setup_samEntry(&readStackST) )
      /*handel memory error*/
   if( setup_samEntry(&refStackST) )
      /*handel memory error*/
   
   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: get first read*/
      if(*readStackST.extraStr == '@')
         break; /*read in first read*/
   } /*Loop: get first read*/

   if(
      mkDepthProfile_edDist(
         &refStackST,
         7,       /*base has q-score of 7 or more*/
         0.7,     /*reads share at least 70% overlap*/
         &distanceStackST,
         &readStackST,/*function uses to get reads*/
         samFILE
      )
   ){
      /*deal with errors*/
   }

   /*DO SOMETHING WITH PROFILE HERE*/

   freeStack_res_edDist(&distanceStackST);
   freeStack_samEntry(&refStackST);
   freeStack_samEntry(&readStackST);

   fclose(samFILE);
   samFILE = 0;

   return 0;
}
```


### read to read edit distance

The read to read edit distance can be found with
  `readCmpDist\_edDist` (fun06). It takes several inputs
  and returns the edit distance between two reads. This
  function is less accurate because it does not attempt
  to re-align the reads. So, differences in the alignments
  can give result in distances greater then the real edit
  distance.

Along with Q-score you can also filter each nucletoide by
  its read depth. However, this requires you to first find
  the read depths with mkDepthProfile\_edDist (see last
  section) or adding each read one by one with
  addReadToDepth\_edDist (see last section).

- Input:
  - samEntry structure pointer with the query read
  - samEntry structure pointer with the reference read to
    compare to
    - must be mapped to the same reference as the query
  - minimum length to count an indel
    - any indel beneath this length is ingored
  - minimum quality score to count a base in the edit
    distance
    - applies to both query and reference read
  - minimum overlap (as percentage [0 to 1]) needed to
    find the edit distance between the two reads
  - minimum depth seen in other reads for this position
    (the depths are stored in the `res\_edDist` structure)
  - size of one window
    - edDist counts distance by windows (ex 500 base
      windows) and total distance
  - `res\_edDist` structure to store the extra results in
    - ex: window results, overlap, and indel counts
    - includes depthAryUI, which is filled by
      mkDepthProfile\_edDist, addReadToDepth\_edDist, or is
     0

- Returns:
  - edit distance (0 or greater) for no errors
  - def\_noSeq\_edDist (is a negative number) if the
    read has no sequence
    - this would be a secondary alignment
  - def\_noMap\_edDist (is a negative number) if the
    read did not map to the reference genome
  - def\_diffRef\_edDist if both query and reference were
    mapped to a different genome (reference)
  - def\_noOverlap\_edDist if query and reference do not
    have enough overlap

Get edit distances for all other reads to the first read
  in a sam file.

```
#ifdef PLAN9
   #include <u.h>
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "edDist.h"

int
main(
){
   signed long distSL = 0;
   
   struct res_edDist distanceStackST;
   struct samEntry readStackST;
   struct samEntry refStackST
   
   FILE *samFILE = fopen("reads.sam", "r");
   
   init_res_edDist(&distanceStackST);
   init_samEntry(&refStackST);
   init_samEntry(&readStackST);
   
   if( setup_samEntry(&readStackST) )
      /*handel memory error*/
   if( setup_samEntry(&refStackST) )
      /*handel memory error*/
   
   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: get first read*/
      if(*readStackST.extraStr == '@')
         break; /*read in first read*/
   } /*Loop: get first read*/

   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: read sam file*/
      distSL =
         readCmpDist_edDist(
            &readStackST,/*read to get distance for*/
            &refStackST, /*reference sequence mapped to*/
            5,         /*indels be at least 5 bases long*/
            7,         /*base has q-score of 7 or more*/
            0.7,      /*reads share at least 70% overlap*/
            0,         /*do not apply read depths*/
            500,       /*500 bases per window*/
            &distanceStackST
         );
   
      if(distSL == def_noSeq_edDist)
         continue; /*header or secondary alignment*/
      else if(distSL < 0)
         continue; /*unmapped read */
   
      printf(
         "distance for %s: %li\n",
         readStackST.qryIdStr,  /*read id*/
         distSL                 /*edit distance*/
      );
   } /*Loop: read sam file*/
   
   freeStack_res_edDist(&distanceStackST);
   freeStack_samEntry(&refStackST);
   freeStack_samEntry(&readStackST);

   fclose(samFILE);
   samFILE = 0;

   return 0;
}
```

# Printing and precentages

The percDist\_edDist (fn08) finds the percent distance
  between two reads. The formula is
  `100 * (distance / (overlap * expected error rate)`.
  Not very great, but from my experince percentages are
  better to compare then edit distances, because edit
  distances get larger the longer your read is. While
  percentages do not.

- Input:
  - `res_edDist` structure pointer with edit distance to
     find percentage for
   - window size used with the edit distence functions
  - expect percent error [0 to 1].
- Output:
  - integer with percentage (rounded to a whole percent)
  - also modifies percentage variables in the input
    `res_edDist ` structure (probMaxWinUI and probUI)

If you want to print out your results like edDist, you
  can use `phead_edDist` (fun11) to print out the header.
  The only input is a FILE pointer to the file to print
  to. The output is the edDist header printed to that
  file.

To print the stats in a `res_edDist` structure  you can
  use `pdist_edDist`.

- Input:
  - `res_edDist` structure pointer to print to file
  - c-string with the read id/name
  - c-string with the reference id/name
  - FILE pointer to file to print to
- Output:
  - prints stats in the `res_edDist` strucure to the
    output file

# Full examples

## reference based edit distance example

A full example for this comparing the reads to their
  mapped reference might look like:

```
#ifdef PLAN9
   #include <u.h>
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "seqST.h"
#include "samEntry.h"
#include "edDist.h"

int
main(
){
   signed long distSL = 0;
   
   struct seqST refStackST
   struct res_edDist distanceStackST;
   struct samEntry readStackST;
   
   FILE *refFILE = fopen("reference.fa", "r");
   FILE *samFILE = fopen("reads.sam", "r");
   
   init_res_edDist(&distanceStackST);
   init_seqST(&refStackST);
   init_samEntry(&readStackST);
   
   if( setup_samEntry(&readStackST) )
      /*handel memory error*/
   
   if( getFa_seqST(refFILE, &refStackST) )
      /*handle the file error or memory error*/
   
   fclose(refFILE);
   refFILE = 0;
   
   phead_edDist(stdout);

   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: read sam file*/
      distSL =
         dist_edDist(
            &readStackST,/*read to get distance for*/
            &refStackST, /*reference sequence mapped to*/
            5,         /*indels be at least 5 bases long*/
            7,         /*base has q-score of 7 or more*/
            500,       /*500 bases per window*/
            &distanceStackST
         );
   
      if(distSL == def_noSeq_edDist)
         continue; /*header or secondary alignment*/
      else if(distSL < 0)
         continue; /*unmapped read */
   
      pdist_edDist(
         readStackST.qryIdStr,
         readStackST.refIdStr,
         stdout
      );
   } /*Loop: read sam file*/
   
   freeStack_res_edDist(&distanceStackST);
   freeStack_seqST(&refStackST);
   freeStack_samEntry(&readStackST);

   fclose(samFILE);
   samFILE = 0;

   return 0;
}
```
## read to read edit distance example (depths ignored)

A full example for this comparing the first read to all
  other reads might look like: This does not use read
  depth in the comparison

```
#ifdef PLAN9
   #include <u.h>
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "edDist.h"

int
main(
){
   signed long distSL = 0;
   
   struct res_edDist distanceStackST;
   struct samEntry readStackST;
   struct samEntry refStackST
   
   FILE *samFILE = fopen("reads.sam", "r");
   
   init_res_edDist(&distanceStackST);
   init_samEntry(&refStackST);
   init_samEntry(&readStackST);
   
   if( setup_samEntry(&readStackST) )
      /*handel memory error*/
   if( setup_samEntry(&refStackST) )
      /*handel memory error*/
   
   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: get first read*/
      if(*readStackST.extraStr == '@')
         break; /*read in first read*/
   } /*Loop: get first read*/


   phead_edDist(stdout);

   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: read sam file*/
      distSL =
         readCmpDist_edDist(
            &readStackST,/*read to get distance for*/
            &refStackST, /*reference sequence mapped to*/
            5,        /*indels be at least 5 bases long*/
            7,        /*base has q-score of 7 or more*/
            0.7,      /*reads share at least 70% overlap*/
            0,        /*ingore read depth*/
            500,      /*500 bases per window*/
            &distanceStackST
         );
   
      if(distSL == def_noSeq_edDist)
         continue; /*header or secondary alignment*/
      else if(distSL < 0)
         continue; /*unmapped read */
   
      pdist_edDist(
         readStackST.qryIdStr,
         refStackST.qryIdStr,
         stdout
      );
   } /*Loop: read sam file*/
   
   freeStack_res_edDist(&distanceStackST);
   freeStack_samEntry(&refStackST);
   freeStack_samEntry(&readStackST);

   fclose(samFILE);
   samFILE = 0;

   return 0;
}
```

## read to read edit distance example (depths used)

A full example for this comparing the first read to all
  other reads might look like. Here the read depth for
  each base is used in the comparison.

```
#ifdef PLAN9
   #include <u.h>
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "samEntry.h"
#include "edDist.h"

int
main(
){
   signed long distSL = 0;
   
   struct res_edDist distanceStackST;
   struct samEntry readStackST;
   struct samEntry refStackST
   
   FILE *samFILE = fopen("reads.sam", "r");
   
   init_res_edDist(&distanceStackST);
   init_samEntry(&refStackST);
   init_samEntry(&readStackST);
   
   if( setup_samEntry(&readStackST) )
      /*handel memory error*/
   if( setup_samEntry(&refStackST) )
      /*handel memory error*/
   
   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: get first read*/
      if(*readStackST.extraStr == '@')
         break; /*read in first read*/
   } /*Loop: get first read*/


    /*get read depths; this is not requried*/
    if(
      mkDepthProfile_edDist(
         &refStackST,
         7,       /*base has q-score of 7 or more*/
         0.7,     /*reads share at least 70% overlap*/
         &distanceStackST,
         &readStackST,/*function uses to get reads*/
         samFILE
      )
   ) /*handle error here*/

   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: move to first read*/
      if(*readStackST.extraStr == '@')
         break; /*read in first read*/
   } /*Loop: move to first read*/

   phead_edDist(stdout);

   while(! get_samEntry(&readStackST, samFILE) )
   { /*Loop: read sam file*/
      distSL =
         readCmpDist_edDist(
            &readStackST,/*read to get distance for*/
            &refStackST, /*reference sequence mapped to*/
            5,        /*indels be at least 5 bases long*/
            7,        /*base has q-score of 7 or more*/
            0.7,      /*reads share at least 70% overlap*/
            10,       /*at least 10 reads support*/
            500,      /*500 bases per window*/
            &distanceStackST
         );
   
      if(distSL == def_noSeq_edDist)
         continue; /*header or secondary alignment*/
      else if(distSL < 0)
         continue; /*unmapped read */
   
      pdist_edDist(
         readStackST.qryIdStr,
         refStackST.qryIdStr,
         stdout
      );
   } /*Loop: read sam file*/
   
   freeStack_res_edDist(&distanceStackST);
   freeStack_samEntry(&refStackST);
   freeStack_samEntry(&readStackST);

   fclose(samFILE);
   samFILE = 0;

   return 0;
}
```
