/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' hirschberg SOF: Start Of File
'   o header:
'     - forward declartions and guards
'   o .c fun01: alnOneNt_hirschberg
'     - align a single base to a sequence
'   o .c fun02: forwardScore_hirschberg
'     - scores the input query and reference sequence in a
'        forward direction
'   o .c fun03: reverseScore_hirschberg
'     - scores input query and reference sequence
'       backwards for a hirschberg
'   o .c fun04: recursive_hirschberg
'     - does the recurisve part of a hirschberg alignment
'   o .c fun05: alnToSam_hirschberg
'     - converts an alignment from recursive_hirschberg to
'       a sam file entry
'   o fun06: aln_hirschberg
'     - does hischberg alignment
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - forward declarations and guards
\-------------------------------------------------------*/

#ifndef HIRSCHBERG_ALGINMENT_H
#define HIRSCHBERG_ALGINMENT_H

struct seqST;
struct samEntry;
struct alnSet;

#define def_memErr_hirschberg ( (signed long) ( ((unsigned long) -1) - (((unsigned long) -1) >> 1) ) )
   /*this should give the high (negative) bit of a signed
   `  long
   */

/*-------------------------------------------------------\
| Fun06: aln_hirschberg
|  - does hischberg alignment
| Input:
|  - refST:
|    o Pointer to seqST structure with the reference
|      sequence and the first base (offsetUL) and last
|      base (endAlnUL) to align (both are index 0)
|  - qryST:
|    o Pointer to seqST structure with the query
|      sequence and the first base (offsetUL) and last
|      base (endAlnUL) to align (both are index 0)
|  - samSTPtr:
|    o pointer to samEntry struct to store alignment in
|  - settings:
|    o Pointer to alnSet structure with the settings for
|      the alignment
| Output:
|  - Returns:
|    o score of alignment for success
|    o def_memErr_hirschberg for memory errors
\-------------------------------------------------------*/
signed long
aln_hirschberg(
  struct seqST *qryST, /*qeury sequence to align*/
  struct seqST *refST, /*reference sequence to align*/
  struct samEntry *samSTPtr,/*gets alignment*/
  struct alnSet *settings   /*settings for alignment*/
);

#endif
