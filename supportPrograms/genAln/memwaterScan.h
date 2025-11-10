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
'   o fun06: qryByRefSort_aln_memwaterScan
'     - sorts the query coordinates in a aln_memwaterScan
'       structure by the reference coordinate
'     - the sorting is first smallest query state, then
'       highest score, then longest (ending base) length
'   o fun07 memwaterScan:
'     - performs a memory efficent Smith Waterman scan
'       (keep best alignment for each query/reference base)
'       alignment on a pair of sequences
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
   `   reference base
   ` Also the best score for each reference base
   */
   signed long *refStartArySL;
   signed long *refEndArySL;
   signed long *refScoreArySL;

   /*alignment coordiantes (has both reference and query
   `   coordiantes, convert to individual coordinates with
   `   seqToIndex) of the best alignment for each
   `   query base
   ` Also the best score for each query base
   */
   signed long *qryStartArySL;
   signed long *qryEndArySL;
   signed long *qryScoreArySL;
   signed int qrySizeSI; /*maxiumum bases in query before
                         `  need to do a resize
                         */
   /*general variables used by memwaterScan and its
   `  supporting functions
   */
   signed long *scoreArySL; /*scoring row*/
   signed long *indexArySL; /*start index of each row*/
   signed char *dirArySC;   /*direction row*/
   signed int rowSizeSI;    /*maximum bytes before need to
                            `  resize scoreArySL,
                            `  indexArySL, dirArySC,
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
|   - qryBl:
|     o 1 if index is for the query array
|     o 0 if index is for the reference array
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
|     o frees alnSTPtr; you must set the pointer to 0/null
|   - Returns:
|     o 0 for success
|     o 1 if no position at index
|     o -1 for out of bounds error
\-------------------------------------------------------*/
signed char
getCoord_aln_memwaterScan(
   signed int baseSI,
   signed char qryBl,         /*1: get from query array*/
   signed int *refStartSIPtr, /*gets reference start*/
   signed int *refEndSIPtr,   /*gets reference end*/
   signed int *qryStartSIPtr, /*gets query start*/
   signed int *qryEndSIPtr,   /*gets query end*/
   struct aln_memwaterScan *alnSTPtr
);

/*-------------------------------------------------------\
| Fun06: qryByRefSort_aln_memwaterScan
|   - sorts the query coordinates in a aln_memwaterScan
|     structure by the reference coordinate
|   - the sorting is first smallest query state, then
|     highest score, then longest (ending base) length
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
qryByRefSort_aln_memwaterScan(
   struct aln_memwaterScan *alnSTPtr
);

/*-------------------------------------------------------\
| Fun07: memwaterScan
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

#endif
