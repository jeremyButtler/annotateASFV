/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' kmerFun SOF: Start Of File
'   - has functions to converts a nucleotide index
'     sequence (seqToIndex_alnSet) into to a kmer sequence
'   o header:
'     - guards
'   o fun01: seqToKmer_kmerFun
'     - convert sequence to a kmer list
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards
\-------------------------------------------------------*/

#ifndef KMER_FUNCTIONS_H
#define KMER_FUNCTIONS_H

/*these are from kmerBit, but are needed when working
`   with kmerFun functions
*/
#define def_bitsPerKmer_kmerFun 2 /*do not change*/
#define def_endKmers_kmerFun -2
#define def_noKmer_kmerFun -1

/*-------------------------------------------------------\
| Fun01: seqToKmer_kmerFun
|   - convert sequence to a kmer list
| Input:
|   - seqStr:
|     o c-string with sequence to convert
|   - lenSeqSI:
|     o length of sequence to convert
|   - kmerArySIPtr:
|     o signed int array pionter to hold kmer sequence
|   - sizeSIPtr:
|     o signed int pointer with current kmerArySIPtr size
|   - kmerLenUC:
|     o length of one kmer (index 1)
| Output:
|   - Modifies:
|     o kmerArySIPtr to have forward kmers (resized if
|       needed)
|     o sizeSIPtr to have resized array length
|   - Returns:
|     o number kmers in sequence
|     o 0 for memory error
\-------------------------------------------------------*/
signed int
seqToKmer_kmerFun(
   signed char *seqStr,        /*sequence to convert*/
   signed int lenSeqSI,        /*sequence length*/
   signed int **kmerArySIPtr,  /*gets kmer sequence*/
   signed int *sizeSIPtr,      /*current array sizes*/
   unsigned char kmerLenUC     /*size of one kmer*/
);

#endif
