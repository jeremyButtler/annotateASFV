# Goal:

Describe how to use the kmerFind library to find the
  locatoin of small (10 to 40 base pair) sequences in
  larger sequences.

- Files: kmerFind.h and kmerFind.c

- Depdendencies:
  - memWater
  - alnSet
  - alnDefs.h
  - seqST from genBio
  - fileFun from genLib
  - kmerBit.h from genBio
  - shellSort from genLib
  - ulCp from genLib
  - genLine.h from genLib
  - genMath.h from genLib (no .c file used)
- Hidden Dependencies:
  - indexToCoord
  - base10str from genLib

You should be familar with the seqST library from genBio
  (see ../genBio/using\_seqST.md).

# Overview:

The kmerFind library is good at searching for the
  positions of small sequences, such as primers in longer
  sequences. However, larger sequences are often very
  inefficent. As a rule, it will be faster, but less
  sensative then a Waterman, but will be able to detect
  smaller alignments (ex 10 to 20 base pair) that minimap2
  would fail to align. Larger alignments become very
  ineffecient, both in memory and time.

This is the library behind the primFind program in
  bioTools.

There are two structures in kmerFind. The first is
  tblST\_kmerFind, which holds the sequence I am aligning
  the small sequences. The second is refST\_kmerFind,
  which holds the small sequences I am aligning.

- Method:
  1. convert longer sequence into kmers
  2. convert short sequences (aligning) into kmers
     - both forward and reverse sequences are used
  3. break longer sequence into windows
     - one window is the size of the longest of the
       short sequences
  4. count kmers shared between longest sequence and
     short sequences for two windows
     - if have 40% or more shared kmers, move to step 5
     - else ignore and move to step 6
  5. if enough matches, align matching short sequence to
     two windows in the longer sequence
     - if high enough score, is a match, record
       - for pratical purposes, I only record the first
         highest score
       - you could modify `fxFindPrims_kmerFind()` to
         record all hits
     - else not match, ignore
  6. move one window and repeat steps 1 to 6 untill no
     more windows

## tblST\_kmerFind structure

The tblST\_kmerFind structure holds the sequence I am
  aliging the smaller seqeunces to. As a rule, you should
  not have to work with this directly.

- Variables:
  1. lenKmerUC: number of nucleotides in one kmer
  2. tblSI: signed int array (really table) with the
     number of times each kmer appears in a window
     - this is a kmer table, so it means it is the size of
       all possible kmers
     - its length/size is stored in lenTblUI
  3. kmerArySI: the kmers in the current windows I am
     comparing (default is two windows) for the longer
     sequence
     - number of stored kmers is in numKmerUI
  4. ntInWinUI: this is the maximum number of kmers I
     can compare at one time
  5. rmNtUI: number of bases to shift by after each
     comparison
  6. kmerMaskUL: mask for clearing old kmers from each
     signed int (do not mess with)
  7. lastKmerUL: last kmer I was working on
     - used in finding the next kmer
     - the number of bases in lastKmerUL is stored in
       lenLastKmerUL
  8. seqPosUL: position at in the longest sequence
  9. seqSTPtr: seqST structure pointer with the longest
     sequence (used in aligning bases)

### tblST\_kmerFind general functions

For using a tblST\_kmerFind structure, the more general
  functions in kmerFind will handel a lot of the details
  for you. So, the main steps are initialize, setup, and
  freeing at the end.

To initialize a tblST\_kmerFind structure, use
  the `init_tblST_kmerFind()` function. The input is a
  tblST\_kmerFind structure pointer to initialize. Only
  call this function once.

The setup step allocated the memory for the tblST\_kmerST
  structure.

- Input:
  1. tblST\_kmerFind structure to initialize
  2. number of nucloetides in a kmer
- Returns:
  1. 0 for no errors
  2. def\_memErr\_kmerFind for memory errors

For freeing, use the `freeStack_tblST_kmerFind()` function
  to free variables inside the structure. You can use
  the `freeHeap_tblST_kmerFind()` function to free a heap
  allocted structure (`freeStack_tblST_kmerFind()` +
  `free()`). Just remember to set your pointer to 0/null
  when finished.

- Input (stack and heap):
  1. tblST\_kmerFind structure pointer to free

**Example**

```
struct tblST_kmerFind tblStackST;
init_tblST_kmerFind(&tblStackST);

if( setup_tblST_kmerFind(&tblStackST) )
   /*deal with memory error*/

/*use your structure*/

freeStack_tblST_kmerFind(&tblStackST);
```

**misc functions**

In addition to the initialization, setup, and freeing
  functions you can also use `blank_tblST_kmerFind()` and
  `qckBlank_tblST_kmerFind()` to blank (prepares for
  another sequence; set used values to 0) a
  tblST\_kmerFind structure. The qckBlank lazy sets all
  lengths to 0, while the blank does a full overwrite.

As a rule you should not need to blank,
  since `fxFindPrims_kmerFind()` will do this for you.

The next step after setup, would be to pepare your
  tblST\_kmerFind structure for alignment. This would use
  the `prep_tblST\_kmerFind()` function, however,
  the `tsvToAry_refST_kmerFind()` or
  the `faToAry_refST_kmerFind()` functions used to get
  your refST\_kmerFind will auto call this for you. So,
  you can skip this step.

After preparation, the next step would be adding in the
  sequence to search. There is no function to do this step
  indvidually, so if are not
  using `fxFindPrims_kmerFind()`, you will have to do
  this manually. The `fxFindPrims_kmerFind()` function
  handels this step for you.

## refST\_kmerFind

The refST\_kmerFind structure holds the sequences to find
  aligment positions for. This would be the short
  sequences. Naming it ref was a bit of a mistake, but
  made sense because I was fining primer or spoligotype
  positions in a read.

- Variables:
  1. lenKmerUC: number of nucleotides in one kmer
  2. minKmersUI: minimum number of kmers needed to do
     an alignment
  3. forKmerArySI: an array with each unique kmer in the
     forward sequence
  4. forRepAryUI: number of times each kmer in the forward
     sequence was repeated
  5. maxForScoreF: maximum alignment score possible for
     the foward sequence
  6. revKmerArySI: an array with each unique kmer in the
     reverse complement sequence
  7. forRepAryUI: number of times each kmer in the reverse
     complement sequence was repeated
  8. maxRevScoreF: maximum alignment score possible for
     the reverse complement sequence
  9. forSeqST: seqST structure with the forward sequence
  10. revSeqST: seqST structure with the reverse
      complement sequence
  11. matSI: the reverse/forward primer this pairs with

There are also several lenght and size variables for the
  arrays.

The kmer arrays (forKmerArySI and revKmerArySI) both end
  with def\endKmers\_kmerBit.

### refST\_kmerFind initialize, setup, and free

As a rule you should not do much to initialize and setup
  this structure, instead you should use
  the `tsvToAry_refST_kmerFind()` or
  the `faToAry_refST_kmerFind()` functions to build your
  refST\_kmerFind structure and prepare your
  tblST\_kmerFind structure. So, I am going to breifly
  describe each indvidual step. Look later for more
  details.

If you plan to use your own functions, you can initialize
  with the `init_refST_kmerFind()` function. You can then
  setup a structure with the `setup_refST_kmerFind()`
  function. You can then use
  the `addSeqToRefST_kmerFind()` function to add a
  sequence to a refST\_kmerFind structure.

Adding sequences will require you to include a setup
  tblST\_kmerFind structure. The return value was designed
  for adding to arrays, were it represents the longest
  sequence you have added. You provide the length of
  the last longest sequence and this decides if the new
  sequence is longer. It then returns the maximum.

After all sequences are added, then you can prepare the
  tblST\_kmerFind structure for scanning with
  the `prep_tblST\_kmerFind()` function.

Finally, you can free the variables in a single
  refST\_kmerFind structure with
  the `freeStack_refST_kmerFind()` function. For a heap
  variable use `freeHeap_refST_kmerFind()`. For heap
  frees, remember to set the pointer to null/0 when done.

- Input: (stack and heap)
  1. pointer to refST\_kmerFind structure to free

You can also free heap allocated refST\_kmerFind arrays
  with the `freeHeapAry_refST_kmerFind()` function. Just
  remember to set your pointer to null/0 when finished.

- Input:
  1. pointer to heap allocated refST\_kmerFind structure
     array to free
  2. length of the array

**Example (using 5mers)**

This example would be for if you are trying to avoid
  the `tsvToAry_refST_kmerFind()` or
  the `faToAry_refST_kmerFind()` functions.

This example has not been tested, but should give you a
  rough idea.

```
/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec01: variable declarations
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

signed int errSI = 0;
struct tblST_kmerFind tblStackST;
struct refST_kmerFind *refHeapAryST = 0;
signed int lenSI = 0;
unsinged int longestSeqUI = 0;

struct seqST inSeqStackST;
struct alnSet alnSettingsStackST;

FILE *inFILE = fopen("primers.fa", "r");

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec02: initialize and setup
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

init_tblST_kmerFind(&tblStackST);
init_seqST(&inSeqStackST);
init_alnSet(&alnSettingsStackST);

if(setup(5, &tblSTStack))
   goto err_main;

refHeapAryST = malloc(5 * sizeof(struct refHeapAryST));
if(! refHeapAryST)
   goto err_main;

if(! inFILE)
   goto err_main;

for(lenSI = 0; lenSI < 5; ++lenSI)
{ /*Loop: initialize and setup*/
   init_refST_kmerFind(&refHeapAryST[lenSI]);
   if( setup_refST_kmerFind(5, &refHeapAryST[lenSI]) )
      goto err_main;
} /*Loop: initialize and setup*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec03: add sequences to search for to the refST_kmerFind
^        structure
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

for(errSI = 0; errSI < lenSI; ++errSI)
{ /*Loop: get five primers from the file*/
   if(getFa_seqST(inFILE, &inSeqStackST))
      goto err_main;

   longestSeqUI =
      addSeqToRefST_kmerFind( 
         &tblStackST,
         &refHeapAryST[errSI], /*add sequence to*/
         &inSeqStackST,        /*sequence to add*/
         0.4,        /*min % of kmers matching .4 = 40%*/
         longestSeqUI,        /*current longest sequence*/
         &alnSettingsStackST
      );

   if(! longestSeqUI)
      goto err_main; /*means memory error*/
} /*Loop: get five primers from the file*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec04: setup the tblST_kmerFind structure
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

if(
   prep_tblST_kmerFind(
      &tblStackST,
      1.0f,    /*1 window is 2x longest sequence length*/
      1.0f,    /*move 1 seq length; 1/2 window per shift*/
      longestSeqUI /*longest sequence length*/
   )
) goto err_main; /*means memory error*/

/*at this point you are ready to call
`  fxFindPrims_kmerFind()
*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec05: deal with errors and clean up
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

errSI = 0;
goto ret_main;

err_main:;
   errSI = 1;
   goto ret_main;

ret_main:;
   freeStack_tblST_kmerFind(&tblStackST);
   freeStack_seqST(&inSeqStackST);
   freeStack_alnSet(&alnSettingsStackST);

   if(inFILE)
      fclose(inFILE);
   inFILE = 0;

   if(refHeapAryST)
      freeHeapAry_refST_kmerFind(refHeapAryST, lenSI);
   refHeapAryST = 0;

   return errSI;
```

# Functions:

- The functions in kmerFind are:
  - Functions used to setup the tblST\_kmerFind structure
    (see structures)
  - Functions to get sequences from file, which return a
    refST\_kmerFind structure. Functions that 
  - Functions to get alignment coordinates
  - Functions to free your structures
  - Functions you do not need to worry about, but are
    needed for the general functions to work
    - I will not cover these in any depth, best to check
      the .c and .h files
