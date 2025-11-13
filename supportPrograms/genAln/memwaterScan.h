/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' memwaterScan SOF: Start Of Functions
'   - has memory efficent Waterman Smith alignment scan
'     (best alignment for each reference/query base)
'   o header:
'     - included libraries
'   o .h st01: aln_memwaterScan
'     - stores the alignment results from a waterman scan
'   o fun01: blank_aln_memwaterScan
'     - blanks (sets non-array values to 0) an
'   o fun02: init_aln_memwaterScan
'     - initializes (sets arrays to 0 + blanks) an
'       aln_memwaterScan structure
'   o fun03: freeStack_aln_memwaterScan
'     - frees the variables in an aln_memwaterScan struct
'   o fun04: freeHeap_aln_memwaterScan
'     - frees an aln_memwaterScan struct
'   o fun05: getCoord_aln_memwaterScan
'     - gets the coordinates for input base index for an
'       aln_memwaterScan structure
'   o fun07: refCoordSort_aln_memwaterScan
'     - sorts best bases alignments by reference start,
'       query start, reference end (closest first), then
'       and query end (closest first)
'     - all negatives (non-alignment found) are pushed to
'       the end
'   o fun08: filter_memwaterScan
'     - removes low scoring alignments and alignments that
'       are nested alignments
'   o fun09: merge_aln_memwaterScan
'     - merges one aln_memwaterScan struct into another
'       aln_memwaterScan struct
'   o fun10 memwaterScan:
'     - performs a memory efficent Smith Waterman scan
'       (keep best alignment for each query/reference base)
'       alignment on a pair of sequences
'   o fun11: simple_memwaterScan
'     - performs a memory efficent Smith Waterman scan
'       (keep best alignment for each query/reference
'       base) alignment on a pair of sequences
'     - simple means no settings or seqST structure used
'       and no match matrix
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - forward declarations and guards
\-------------------------------------------------------*/

#ifndef MEMORY_EFFICENT_WATERMAN_MATRIX_SCAN_H
#define MEMORY_EFFICENT_WATERMAN_MATRIX_SCAN_H

struct seqST;
struct alnSet;

/*-------------------------------------------------------\
| ST01: aln_memwaterScan
|   - stores the alignment results from a waterman scan
\-------------------------------------------------------*/
typedef struct aln_memwaterScan
{
   /*length of aligned sections of the sequences;
   `  ref=reference, qyr=query
   */
   signed int refLenSI;
   signed int qryLenSI;

   /*offset (starting alignment position) of sequences*/
   signed int refOffsetSI;
   signed int qryOffsetSI;

   /*alignment coordiantes (has both reference and query
   `   coordiantes, convert to individual coordinates with
   `   seqToIndex) of the best alignment for each
   `   reference (0 to refLenSI - 1) and query
   `   (refLenSI to refLenSI + qryLen - 1)
   `   base
   ` Also the best scores for each bases alignment
   */
   signed long *startArySL;
   signed long *endArySL;
   signed long *scoreArySL;
   signed long outSizeSL;
   signed long outLenSL;

   /*general variables used by memwaterScan and its
   `  supporting functions
   */
   signed long *scoreRowSL; /*scoring row*/
   signed long *indexRowSL; /*start index of each row*/
   signed char *dirRowSC;   /*direction row*/
   signed int rowSizeSI;    /*maximum bytes before need to
                            `  resize scoreRowSL,
                            `  indexRowSL, dirRowSC,
                            `  refStartArySL, refEndArySL,
                            `  and refScoreArySL
                            */
}aln_memwaterScan;

/*-------------------------------------------------------\
| Fun01: blank_aln_memwaterScan
|   - blanks (sets non-array values to 0) an
|     aln_memwaterScan struct
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer to blank
| Output:
|   - Modifies:
|     o sets refLenSI, qryLenSI, refOffsetSI, qryOffsetSI
|       in alnsTPtr to 0
\-------------------------------------------------------*/
void
blank_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
);

/*-------------------------------------------------------\
| Fun02: init_aln_memwaterScan
|   - initializes (sets arrays to 0 + blanks) an
|     aln_memwaterScan structure
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer to initialize
| Output:
|   - Modifies:
|     o sets all values to 0
\-------------------------------------------------------*/
void
init_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
);

/*-------------------------------------------------------\
| Fun03: freeStack_aln_memwaterScan
|   - frees the variables in an aln_memwaterScan struct
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer with variable to
|       free
| Output:
|   - Modifies:
|     o frees all arrays and initializes
\-------------------------------------------------------*/
void
freeStack_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
);

/*-------------------------------------------------------\
| Fun04: freeHeap_aln_memwaterScan
|   - frees an aln_memwaterScan struct
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer to free
| Output:
|   - Modifies:
|     o frees alnSTPtr; you must set the pointer to 0/null
\-------------------------------------------------------*/
void
freeHeap_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
);

/*-------------------------------------------------------\
| Fun05: getCoord_aln_memwaterScan
|   - gets the coordinates for input base index for an
|     aln_memwaterScan structure
| Input:
|   - baseSI:
|     o base to get coordinates for
|   - refStartSIPtr:
|     o signed int pointer that gets the reference start
|   - refEndSIPtr:
|     o signed int pointer that gets the reference end
|   - qryStartSIPtr:
|     o signed int pointer that gets the query start
|   - qryEndSIPtr:
|     o signed int pointer that gets the query end
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer with coordinates
| Output:
|   - Modifies:
|     o refStartSIPtr to have reference starting position
|     o refEndSIPtr to have reference ending position
|     o qryStartSIPtr to have query starting position
|     o qryEndSIPtr to have query ending position
|   - Returns:
|     o score for success
|     o 0 if no position at index
|     o -1 for out of bounds error
\-------------------------------------------------------*/
signed char
getCoord_aln_memwaterScan(
   signed long baseSL,        /*index to grab*/
   signed int *refStartSIPtr, /*gets reference start*/
   signed int *refEndSIPtr,   /*gets reference end*/
   signed int *qryStartSIPtr, /*gets query start*/
   signed int *qryEndSIPtr,   /*gets query end*/
   struct aln_memwaterScan *alnSTPtr
);

/*-------------------------------------------------------\
| Fun06: swap_memwaterScan
|   - swaps position in reference array
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pointer to swap positions
|   - firstSI:
|     o first index to swap
|   - secSI:
|     o second index to swap
| Output:
|   - Modifies:
|     o scoreArySL, startArySL, endArySL in alnSTPtr to
|       have firstSI and secSI index's swapped
\-------------------------------------------------------*/
void
swap_memwaterScan(
   struct aln_memwaterScan *alnSTPtr,
   signed int firstSI,
   signed int secSI
);

/*-------------------------------------------------------\
| Fun07: refCoordSort_aln_memwaterScan
|   - sorts best bases alignments by reference start,
|     query start, reference end (closest first), then
|     and query end (closest first)
|   - all negatives (non-alignment found) are pushed to
|     the end
| Input:
|   - alnSTPtr:
|     o aln_memwaterScan struct pionter to sort
| Output:
|   - Modifies:
|     o qryStartArySL in alnSTPtr to be sorted by
|       reference coodrinate
|     o qryEndArySL and qryScoreArySL in alnSTPtr to be
|       kept in sync with qryStartArySL
\-------------------------------------------------------*/
void
refCoordSort_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
);

/*-------------------------------------------------------\
| Fun08: filter_memwaterScan
|   - removes low scoring alignments and alignments that
|     are nested alignments
| Input:
|   - alnSTPtr:
|     o pointer to aln_memwaterScan struct with alignments
|       to filter
|   - minScoreSL:
|     o minimum score for an alignment
|   - minScoreTwoSL:
|     o second minimum score for an alignment
|       * use when you are putting a percent minimum score
|         for minScoreSL
|       * allows for a hard cutoff and a percent cuttoff
| Output:
|   - Modifies:
|     o startArySL, endArySL, and scoreArySL in alnSTPtr
|       to only have high scoring, non-nested alignments
|     o outLenSL to have the new number of alignments
|   - Returns:
|     o number of kept alignments
\-------------------------------------------------------*/
signed long
filter_memwaterScan(
   struct aln_memwaterScan *alnSTPtr,
   signed long minScoreSL,
   signed long minScoreTwoSL
);

/*-------------------------------------------------------\
| Fun09: merge_aln_memwaterScan
|   - merges one aln_memwaterScan struct into another
|     aln_memwaterScan struct
| Input:
|   - mergeSTPtr:
|     o aln_memwaterScan struct to merge into the final
|       structure
|   - outSTPtr:
|     o aln_memwaterScan struct to merge mergeSTPtr into
|   - minScoreSL:
|     o minumum score to keep score in mergeSTPtr
| Output:
|   - Modifies:
|     o variables in outSTPtr to include variables in
|       mergeSTPtr and outSTPtr
|   - Returns:
|     o 0 for success
|     o -1 if reference lengths differ (no way to merge)
|     o -2 for memory errors
\-------------------------------------------------------*/
signed char
merge_aln_memwaterScan(
   struct aln_memwaterScan *mergeSTPtr,
   struct aln_memwaterScan *outSTPtr,
   signed long minScoreSL
);

/*-------------------------------------------------------\
| Fun10: memwaterScan
|   - performs a memory efficent Smith Waterman scan
|     (keep best alignment for each query/reference base)
|     alignment on a pair of sequences
| Input;
|   - qrySeqSTPtr:
|     o pointer to an seqST with the query sequence and
|       index 0 coordinates to start (offsetSL)/end
|       (endAlnSL) the alignment.
|   - refSeqSTPtr:
|     o pointer to an seqST with the reference sequence
|       and index 0 coordinates to start (offsetSL)/end
|       (endAlnSL) the alignment.
|   - alnSTPtr:
|     o pointer to aln_memwaterScan structure to hold the
|       results of the alignment
|   - settings:
|     o pointer to an alnSet structure with the gap open,
|       gap extend, and scoring matrix for the alingment
| Output:
|  - Modifies:
|    o variables in alnSTPtr to have the new alignment
|  - Returns:
|    o score for aligment
|    o negative number for memory errors
\-------------------------------------------------------*/
signed long
memwaterScan(
   struct seqST *qrySTPtr, /*query sequence and data*/
   struct seqST *refSTPtr, /*ref sequence and data*/
   struct aln_memwaterScan *alnSTPtr,/*gets alignment*/
   struct alnSet *settings
);

/*-------------------------------------------------------\
| Fun11: simple_memwaterScan
|   - performs a memory efficent Smith Waterman scan
|     (keep best alignment for each query/reference base)
|     alignment on a pair of sequences
|   - simple means no settings or seqST structure used and
|     no match matrix
| Input;
|   - qrySeqStr:
|     o c-string with query sequence
|   - qryLenSI:
|     o length of qrySeqStr (index 1)
|   - refSeqStr:
|     o c-string with reference sequence
|   - refLenSI:
|     o length of refSeqStr (index 1)
|   - alnSTPtr:
|     o pointer to aln_memwaterScan structure to hold the
|       results of the alignment
| Output:
|  - Modifies:
|    o variables in alnSTPtr to have the new alignment
|  - Returns:
|    o score for aligment
|    o negative number for memory errors
\-------------------------------------------------------*/
signed long
simple_memwaterScan(
   signed char *qrySeqStr,   /*query sequence*/
   signed int lenQrySI,      /*length of query sequence*/
   signed char *refSeqStr,   /*reference sequence*/
   signed int lenRefSI,      /*reference sequence length*/
   struct aln_memwaterScan *alnSTPtr /*gets alignment*/
);

#endif
