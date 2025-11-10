/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' mapRead SOF: Start Of File
'   - maps a read to a reference
'   o header:
'     - defined variables and guards
'   * .h st01: chains_mapRead
'     - holds kmer chains build from kmerChain_mapRead
'      o fun01: blank_chains_mapRead
'        - blanks a chains_mapRead struct
'      o fun02: init_chains_mapRead
'        - initialize a chains_mapRead struct
'      o fun03: freeStack_chains_mapRead
'        - frees all variables in a chains_mapRead struct
'      o fun04: freeHeap_chains_mapRead
'        - frees a chains_mapRead struct
'      o fun05: setup_chains_mapRead
'        - sets up memory for a chains_mapRead struct
'      o fun06: mk_chains_mapRead
'        - makes a heap allocated chains_mapRead struct
'      o fun07: realloc_chains_mapRead
'        - reallocates memroy in a chains_mapRead struct
'   * .h st02: set_mapRead
'     - settings for a aln_mapRead
'     o fun08: blank_set_mapRead
'       - blank a set_mapRead struct
'     o fun09: init_set_mapRead
'       - initializes a set_mapRead struct
'     o fun10: freeStack_set_mapRead
'       - frees variables in a set_mapRead struct
'     o fun11: freeHeap_set_mapRead
'       - frees a set_mapRead struct
'     o fun12: setup_set_mapRead
'       - allocates memory for structs in a set_mapRead
'   * .h st03: ref_mapRead
'     - holds reference variables used in read mapping
'     o fun13: blank_ref_mapRead
'       - blanks (sets to defaults) a ref_mapRead struct
'     o fun14: init_ref_mapRead
'       - initialize a ref_mapRead structure
'     o fun15: freeStack_ref_mapRead
'       - frees variables in a ref_mapRead struct
'     o fun16: freeHeap_ref_mapRead
'       - frees a ref_mapRead struct
'     o fun17: addRef_ref_mapRead
'       - adds a reference sequence to ref_map read struct
'     o fun18: mk_ref_mapRead
'       - makes and initializes a ref_mapRead structure
'     o fun19: getRef_ref_mapRead
'       - gets and stores next reference in ref_mapRead
'     o fun20: getRefBin_ref_mapRead
'       - gets a binary file with reference
'     o fun21: writeRefBin_ref_mapRead
'       - writes a ref_mapRead struct to a binary file
'   * .h st04: aln_mapRead
'     o fun22: blank_aln_mapRead
'       - blanks an aln_mapRead struct
'     o fun23: init_aln_mapRead
'       - initializes an aln_mapRead struct
'     o fun24: freeStack_aln_mapRead
'       - frees variables in an aln_mapRead struct
'     o fun25: freeHeap_aln_mapRead
'       - frees an aln_mapRead struct
'     o fun26: setup_aln_mapRead
'       - setups (memory allocate) an aln_mapRead struct
'     o fun27: mk_aln_mapRead
'       - makes and returns heap allocated aln_mapRead
'   * .c st05: hit_mapRead
'     - private structure to create a linked list of chain
'       index's that I have already aligned
'   o fun28: sortKmerIndex_mapRead
'     - sorts kmer index array by a kmer array
'   o .c fun29: uiFind_mapRead
'     - finds first signed int in a range
'   o .c fun30: findChain_chains_mapRead
'     - finds kmer in range from kmerChain_mapRead return
'   o fun31: kmerChain_mapRead
'     - merge kmers into chains (longer kmers)
'   o .c fun32: numQryChains_mapRead
'     - gets the number of chains assigned to a single
'       query kmer
'   o fun33: mergeChains_mapRead
'     - merges kmer chains from fun31 (kmerChain) into
'       the chain the maximizes base count
'   o fun34: mergeToSam_mapRead
'     - convert output from mergeChains_map read to sam
'       file entry (mapq is set to 0)
'   o fun35: scoreSubAln_mapRead
'     - scores the aligment in samEntry and finds best
'       sub-aligment
'   o fun37: align_mapRead
'     - maps read to reference
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header
|   - defined variables and guards
\-------------------------------------------------------*/

#ifndef MAP_READ_H
#define MAP_READ_H

struct FILE;
struct seqST;
struct samEntry;
struct alnSet;

/*keep these negative, merge Chains returns postive/0 on
`  success
*/
#define def_memErr_mapRead -1
#define def_noChains_mapRead -2
#define def_noScore_mapRead -4
#define def_noAln_mapRead -8
#define def_noStruct_mapRead -16
#define def_badFile_mapRead -32
#define def_EOF_mapRead -64

/*-------------------------------------------------------\
| ST01: chains_mapRead
|   - holds kmer chains build from kmerChain_mapRead
\-------------------------------------------------------*/
typedef struct chains_mapRead
{
   /*chain arrays (start coords, length, and score)*/
   signed int *refArySI;   /*ref kmer starting coords*/
   signed int *qryArySI;   /*query kmer start coords*/
   signed int *chainArySI; /*length of each chain*/
   signed int *scoreArySI; /*maximal bases in a chain*/
   signed int *gapScoreArySI; /*score with gaps*/
      /*has best scores foun in previous round for
      `  a chain when merged with other chain; ant method?
      */
   signed int *nextArySI;
      /*breadcrum trail used to get best score*/

   /*20bytes per chain*/

   /*limiting to ints, because if I really have
   `    over 2 billion items I am in trouble.
   `    At that point I need at least 40Gb of ram and I
   `    excedes the largest solved traveling salesman
   `    problem by many folds (< 5 million cities)
   */
   signed int lenChainSI;  /*length of chain arrays*/
   signed int sizeChainSI; /*size of chain arrays*/

   /*these allow quick jumping/indexing*/
   signed int *indexArySI;/*refAryUI ref kmers index's*/
   signed int lenIndexSI; /*length match/index arrays*/
   signed int sizeIndexSI;/*size of match/index arrays*/
   signed int lastIndexSI;/*last query nt to start chain*/
}chains_mapRead;

/*-------------------------------------------------------\
| ST02: set_mapRead
|   - settings for a aln_mapRead
\-------------------------------------------------------*/
typedef struct set_mapRead
{
   /*determined by reference length*/
   signed int lenArySI[16];  /*length to swith kmer size*/
   unsigned char kmerAryUC[16]; /*kmer size per length*/
   signed int chainLenArySI[16]; /*minimum chain length
                                 `  to keep a chain
                                 */
   unsigned char lenKmersUC;    /*number kmer lengths*/

   /*chain merging*/
   signed int gapSI;   /*gap score for chain merging*/
   signed int matchSI; /*match score for chain merging*/

   /*for gap filling*/
   struct alnSet *alnSetST;
   signed char alnEndsBl; /*align softmasked ends*/

   /*for scoreSub_mapRead*/
   signed char subBl;     /*1: do sub-aligment score*/
   float minScoreF;       /*min % score to keep*/
   float minMatchF;       /*min % matches to keep read*/
   float minPercLenF;     /*min % score to keep*/
   float chainMinNtF;    /*min % bases in chainto keep*/
   float maxLenPercF;     /*max % length between chains*/
}set_mapRead;

/*-------------------------------------------------------\
| ST03: ref_mapRead
|   - holds reference variables used in read mapping
\-------------------------------------------------------*/
typedef struct ref_mapRead
{
   struct seqST *seqSTPtr;      /*reference seq + id*/
   /*cost is ~ 1x, or ~ 2x if have q-score entry*/

   /*settings for mapping*/
   unsigned char lenKmerUC;     /*length of one kmer*/
   signed int minChainLenSI;    /*minimum chain length*/
     /*in most; if not all cases this will be lenKmerUC*/

   /*reference array cost: bytes = 8x reference length*/
   signed int *kmerArySI;     /*reference as kmer seq*/
   signed int sizeKmerSI;     /*size of kmer array*/
   signed int *indexArySI;
      /*convert refSortKmerArySI index to refKmerArySI
      `  index
      */
   signed int sizeIndexSI;    /*sized of sorted arrays*/
   signed int lenSI;          /*number kmers in arrays*/
}ref_mapRead;

/*-------------------------------------------------------\
| ST04: aln_mapRead
|   - holds non-reference variables used in read mapping
\-------------------------------------------------------*/
typedef struct aln_mapRead
{
   signed int *kmerArySI;  /*query forward kmer seq*/
   signed int lenSI;       /*length of arrays*/
   signed int sizeSI;      /*size of arrays*/

   struct dirMatrix *matrixSTPtr;/*matrix for alignment*/
   struct chains_mapRead *chainsSTPtr;
      /*has chains to use at all steps*/
}aln_mapRead;

/*-------------------------------------------------------\
| Fun01: blank_chains_mapRead
|   - blanks a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct to blank
| Output:
|   - Modifies:
|     o refKmerAryUI and qryKmerAryUI in chainsSTPtr to
|       have all values set to -1
|     o chainLenArySI and scoreAryF in chainsSTPtr to have 
|       all values set to 0
|     o indexArySL in chainsSTPtr to be filled with -1's
|     o lenChainSL and lenIndexSL in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
blank_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to blank*/
);

/*-------------------------------------------------------\
| Fun02: init_chains_mapRead
|   - initialize a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead struct to initialize
| Output:
|   - Modifies:
|     o refKmerAryUI, qryKmerAryUI, chainLenArySI,
|       scoreAryF, and indexArySL to be 0 (null)
|     o lenChainSL, sizeChainSL, lenIndexSL, and
|       sizeIndexSL in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
init_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to initialize*/
);

/*-------------------------------------------------------\
| Fun03: freeStack_chains_mapRead
|   - frees all variables in a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead with variables to free
| Output:
|   - Modifies:
|     o frees refKmerAryUI, qryKmerAryUI, chainLenArySI,
|       scoreAryF, and indexArySL and sets to 0 (null)
|     o lenChainSL, sizeChainSL, lenIndexSL, and
|       sizeIndexSL in chainsSTPtr to be 0
\-------------------------------------------------------*/
void
freeStack_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to free vars*/
);

/*-------------------------------------------------------\
| Fun04: freeHeap_chains_mapRead
|   - frees a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to free
| Output:
|   - Modifies:
|     o frees chainsSTPtr (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_chains_mapRead(
   struct chains_mapRead *chainsSTPtr /*to free*/
);

/*-------------------------------------------------------\
| Fun05: setup_chains_mapRead
|   - sets up memory for a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to setup memory
|   - lenChainsSL:
|     o length of chains arrays
|   - lenRefUI:
|     o length of reference
| Output:
|   - Modifies:
|     o allocates lenChainsSL items of memory to
|       refCkmerArySI, qryKmerArySI, chainLenArySI,
|       and scoreAryF in chainsSTPtr
|     o allocates lenRefUI memory to indexArySL in
|       chainsSTPtr
|     o sizeChainSL in chainsSTPtr to be lenChainSL 
|     o sizeIndexSL in chainsSTPtr to be lenIndexSL 
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*to setup*/
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
);

/*-------------------------------------------------------\
| Fun06: mk_chains_mapRead
|   - makes a heap allocated chains_mapRead struct
| Input:
|   - lenChainsSL:
|     o length of chains arrays
|   - lenRefUI:
|     o length of reference
| Output:
|   - Returns:
|     o pointer to initialized and setup chains_mapRead
|       struct
|     o 0 for memory errors
\-------------------------------------------------------*/
struct chains_mapRead *
mk_chains_mapRead(
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
);

/*-------------------------------------------------------\
| Fun07: realloc_chains_mapRead
|   - reallocates memroy in a chains_mapRead struct
| Input:
|   - chainsSTPtr:
|     o pointer to a chains_mapRead to rellocate memory
|   - lenChainsSL:
|     o new length of chains arrays (0 for do not resize)
|   - lenRefUI:
|     o new length of reference (0 for do not resize)
| Output:
|   - Modifies:
|     o refKmerArySI, qryKmerArySI, chainLenArySI,
|       scoreAryF in chainsSTPtr to be reallocated to
|       lenChainSL if lenChainSL >chainsSTPtr->sizeChainSL
|     o indexArySL chainsSTPtr to be reallocated to
|       lenRefUI if lenRefUI  > chainsSTPtr->sizeIndexSL
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
realloc_chains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*reallocate*/
   signed long lenChainsSL,         /*chain array length*/
   unsigned int lenRefUI            /*reference length*/
);

/*-------------------------------------------------------\
| Fun08: blank_set_mapRead
|   - blank a set_mapRead struct
| Input:
|   - setSTPtr:
|     o pointer to set_mapRead struct to blank
| Output:
|   - Modifies:
|     o blanks structs in setSTPtr
\-------------------------------------------------------*/
void
blank_set_mapRead(
   struct set_mapRead *setSTPtr
);

/*-------------------------------------------------------\
| Fun09: init_set_mapRead
|   - initializes a set_mapRead struct
| Input:
|   - setSTPtr:
|     o pointer to set_mapRead struct to initializes
| Output:
|   - Modifies:
|     o sets variables to defaults and pointers to null
\-------------------------------------------------------*/
void
init_set_mapRead(
   struct set_mapRead *setSTPtr
);

/*-------------------------------------------------------\
| Fun10: freeStack_set_mapRead
|   - frees variables in a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer with variables to free
| Output:
|   - Modifies:
|     o frees all pointers and intitializes; defaults/null
\-------------------------------------------------------*/
void
freeStack_set_mapRead(
   struct set_mapRead *setSTPtr
);

/*-------------------------------------------------------\
| Fun11: freeHeap_set_mapRead
|   - frees a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o frees setSTPtr (you must set to null/0)
\-------------------------------------------------------*/
void
freeHeap_set_mapRead(
   struct set_mapRead *setSTPtr
);

/*-------------------------------------------------------\
| Fun12: setup_set_mapRead
|   - allocates memory for structs in a set_mapRead struct
| Input:
|   - setSTPtr:
|     o set_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o alnSetST to be on heap
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_set_mapRead(
   struct set_mapRead *setSTPtr
);

/*-------------------------------------------------------\
| Fun13: blank_ref_mapRead
|   - blanks (sets to defaults) a ref_mapRead structure
|   - use addRef_ref_mapRead to change references
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to blank
| Output:
|   - Modifies:
|     o sets all values in refSTPtr to defaults or for
|       arrays null
\-------------------------------------------------------*/
void
blank_ref_mapRead(
   struct ref_mapRead *refSTPtr
);

/*-------------------------------------------------------\
| Fun14: init_ref_mapRead
|   - initialize a ref_mapRead structure
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to initialize
| Output:
|   - Modifies:
|     o sets pointers to 0/null & other values to default
\-------------------------------------------------------*/
void
init_ref_mapRead(
   struct ref_mapRead *refSTPtr
);

/*-------------------------------------------------------\
| Fun15: freeStack_ref_mapRead
|   - frees variables in a ref_mapRead struct
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer with variables to free
| Output:
|   - Modifies:
|     o frees all heap allocated variabes in refSTPtr
|     o all variabes in refSTPtr to be 0/null (free) or
|       defaults
\-------------------------------------------------------*/
void
freeStack_ref_mapRead(
   struct ref_mapRead *refSTPtr
);

/*-------------------------------------------------------\
| Fun16: freeHeap_ref_mapRead
|   - frees a ref_mapRead struct
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o frees refSTPtr (you must set to 0/null)
\-------------------------------------------------------*/
void
freeHeap_ref_mapRead(
   struct ref_mapRead *refSTPtr
);

/*-------------------------------------------------------\
| Fun17: addRef_ref_mapRead
|   - adds a reference sequence to a ref_map read struct
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to add reference to
|   - seqSTPtr:
|     o seqST struct pointer with reference sequence
|     o must be heap allocated; otherwise messes up
|       freeStack_ref_mapRead
|     o on success, you should set this to null, but do
|       not free
|     o on failure, free seqSTPtr; only assigned on succes
|   - setSTPtr:
|     o set_mapRead struct with settings for alignment
| Output:
|   - Modifies:
|     o seqSTPtr in refSTPtr to point to seqSTPtr
|     o kmerArySI to have refernce kmers
|     o sortKmerArySI to have refernce kmers and be sorted
|     o indexArySI to have index's to convert
|       sortKmerArySI index to kmerArySI index
|     o lenKmerUC and minChainLenSI to bet set to their
|       correct values in setSTPtr
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|     o def_noStruct_mapRead if any input is 0/null
\-------------------------------------------------------*/
signed char
addRef_ref_mapRead(
   struct ref_mapRead *refSTPtr, /*add reference to*/
   struct seqST *seqSTPtr,       /*reference to add*/
   struct set_mapRead *setSTPtr  /*has settings to add*/
);

/*-------------------------------------------------------\
| Fun18: mk_ref_mapRead
|   - makes and initializes a ref_mapRead structure
| Input:
|   - seqSTPtr:
|     o seqST struct pointer with reference sequence
|     o must be heap allocated; otherwise messes up
|       freeStack_ref_mapRead
|     o on success, you should set this to null, but do
|       not free
|     o on failure, free seqSTPtr; only assigned on succes
|   - setSTPtr:
|     o set_mapRead struct with settings for alignment
| Output:
|   - Returns:
|     o 0 for errors
|     o pointer to setup ref_mapRead structure
\-------------------------------------------------------*/
struct ref_mapRead *
mk_ref_mapRead(
   struct seqST *seqSTPtr,       /*reference to add*/
   struct set_mapRead *setSTPtr  /*has settings to add*/
);

/*-------------------------------------------------------\
| Fun19: getRef_ref_mapRead
|   - gets and stores next reference in ref_mapRead struct
| Input:
|   - refSTPtr:
|     o ref_mapRead struct pointer to store reference in
|   - setSTPtr:
|     o set_mapRead struct with settings for alignment
|   - fxBl:
|     o tells if file is fastq (1) or regular fasta (0)
|   - gzBl:
|     o tells if file is gziped (1) or regular fasta (0)
|     o TODO: add gzip support
|   - refFILE:
|     o FILE pointer to file with reference to store
| Output:
|   - Modifies:
|     o seqSTPtr in refSTPtr to have reference sequence
|     o kmerArySI in refSTPtr to have refernce kmers
|     o sortKmerArySI to have refernce kmers and be sorted
|     o indexArySI to have index's to convert
|       sortKmerArySI index to kmerArySI index
|     o lenKmerUC and minChainLenSI to bet set to their
|       correct values in setSTPtr
|     o refFILE to be on next entry in fastx file
|   - Returns:
|     o 0 for errors
|     o def_EOF_mapRead if hit end of file
|     o def_memErr_mapRead for memory errors
|     o def_noStruct_mapRead if setSTPtr or refSTPtr is 0
|     o def_badFile_mapRead for invalid fastx entries
|       * also for gz files (currently not supported)
\-------------------------------------------------------*/
signed char
getRef_ref_mapRead(
   struct ref_mapRead *refSTPtr,/*gets reference seq*/
   struct set_mapRead *setSTPtr,/*has settings to add*/
   signed char fxBl,            /*1: fastq; 0: fasta*/
   signed char gzBl,            /*1: gzipped; 0: not*/
   void *refFILE                /*fx file with reference*/
);

/*-------------------------------------------------------\
| Fun20: getRefBin_ref_mapRead
|   - gets a binary file with reference
| Input:
|   - refSTPtr:
|     o ref_mapRead struct to add reference to
|   - refFILE:
|     o FILE pointer to binary file with reference
|       sequence, kmers, sorted index's, minimum chain
|       length, and kmer length
|     o Format:
|       1: kmer length [of one kmer] (unsigned char)
|       2: minimum kmers in a chain (signed int)
|       3: length of reference id (name) (signed int)
|       4: reference id (signed char)
|       5: length of reference sequence (signed int)
|       6: reference sequence (signed char)
|       7: length of kmer sequence (signed int)
|       8: reference kmer sequence (signed int)
|       9: sorted index's for kmer sequence (signed int)
| Output:
|   - Modifies:
|     o refSTPtr to have values in refFILE
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|     o def_noStruct_mapRead if no structure input
|     o def_badRead_mapRead if no file input or if had
|       file error
\-------------------------------------------------------*/
signed char
getRefBin_ref_mapRead(
   struct ref_mapRead *refSTPtr,/*gets reference*/
   void *refFILE               /*bin file with reference*/
);

/*-------------------------------------------------------\
| Fun21: writeRefBin_ref_mapRead
|   - writes a ref_mapRead struct to a binary file
| Input:
|   - refSTPtr:
|     o ref_mapRead struct to write
|   - refFILE:
|     o FILE pointer to write contents of refSTPtr to
| Output:
|   - Prints:
|     o refSTPtr to file
|       1: kmer length [of one kmer] (unsigned char)
|       2: minimum kmers in a chain (signed int)
|       3: length of reference id (name) (signed int)
|       4: reference id (signed char)
|       5: length of reference sequence (signed int)
|       6: reference sequence (signed char)
|       7: length of kmer sequence (signed int)
|       8: reference kmer sequence (signed int)
|       9: sorted index's for kmer sequence (signed int)
|   - Returns:
|     o 0 for no errors
|     o def_noStruct_mapRead if refSTPtr is 0/null
|     o def_badFile_mapRead if refFILE is 0/null
\-------------------------------------------------------*/
signed char
writeRefBin_ref_mapRead(
   struct ref_mapRead *refSTPtr,/*reference to write*/
   void *refFILE              /*file to save refSTPtr to*/
);

/*-------------------------------------------------------\
| Fun22: blank_aln_mapRead
|   - blanks an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer to blank
| Output:
|   - Modifies:
|     o first element kmerArySI in alnSTPtr to be
|       def_endKmers_kmerBit 
|     o lenSI in alnSTPtr to be 0
|     o dirMatrix in alnSTPtr with blank_dirMatrix
|     o chainsSTPtr in alnSTPtr with blank_chains_mapRead
\-------------------------------------------------------*/
void
blank_aln_mapRead(
   struct aln_mapRead *alnSTPtr
);

/*-------------------------------------------------------\
| Fun23: init_aln_mapRead
|   - initializes an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer to initialize
| Output:
|   - Modifies:
|     o kmerArySI, matrixSTPtr, and chainsSTPtr in
|       alnSTPtr to be 0/null
|     o lenSI and sizeSI in alnSTPtr to be 0
\-------------------------------------------------------*/
void
init_aln_mapRead(
   struct aln_mapRead *alnSTPtr
);

/*-------------------------------------------------------\
| Fun24: freeStack_aln_mapRead
|   - frees variables in an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer with variables to free
| Output:
|   - Modifies:
|     o frees kmerArySI, matrixSTPtr, and chainsSTPtr
|       in alnSTPtr
|     o sizeSI and lenSI in alnSTPtr to be 0
\-------------------------------------------------------*/
void
freeStack_aln_mapRead(
   struct aln_mapRead *alnSTPtr
);

/*-------------------------------------------------------\
| Fun25: freeHeap_aln_mapRead
|   - frees an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer to free
| Output:
|   - Modifies:
|     o frees alnSTPtr (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_aln_mapRead(
   struct aln_mapRead *alnSTPtr
);

/*-------------------------------------------------------\
| Fun26: setup_aln_mapRead
|   - setups (memory allocate) an aln_mapRead struct
| Input:
|   - alnSTPtr:
|     o aln_mapRead struct pointer to setup
| Output:
|   - Modifies:
|     o allocates memory for kmerArySI, matrixSTPtr, and
|       chainsSTPtr in alnSTPtr
|     o sizeSI in alnSTPtr to be size of memory
|     o lenSI in alnSTPtr to be 0
|   - Returns:
|     o 0 for no errors (or no structure input)
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
setup_aln_mapRead(
   struct aln_mapRead *alnSTPtr
);

/*-------------------------------------------------------\
| Fun27: mk_aln_mapRead
|   - makes and returns heap allocated aln_mapRead struct
| Input:
| Output:
|   - Returns:
|     o pointer to new initialized and setup aln_mapRead
|       struct
|     o 0 for memory errors
\-------------------------------------------------------*/
struct aln_mapRead *
mk_aln_mapRead(
   void
);

/*-------------------------------------------------------\
| Fun28: sortKmerIndex_mapRead
|   - sorts kmer index array by a kmer array
| Input:
|   - refKmerArySI:
|     o signed int array of reference kmers to use for
|       sorting
|   - indexArySI:
|     o signed int array to have sorted kmer index's
|     o needs to be same length as lenSI
|   - lenSI:
|     o length of kmer array (index 1)
| Output:
|   - Modifies:
|     o indexArySI to have index's sorted by refKmerArySI
|       - also sets index's for indexArySI
\-------------------------------------------------------*/
void
sortKmerIndex_mapRead(
   signed int *kmerArySI, /*kmers to make sorted index's*/
   signed int *indexArySI,/*get kmerArySI sorted index's*/
   signed int lenSI       /*array length (index 1)*/
);

/*-------------------------------------------------------\
| Fun31: kmerChain_mapRead
|   - merge kmers into chains (longer kmers)
| Input:
|   - chainsSTPtr:
|     o pointer to chainsSTPtr struct to add kmer chains
|       to
|   - minChainLenSI:
|     o minimum number of kmers to keep a chain
|   - kmerUC:
|     o length of one kmer
|   - qryKmerArySI:
|     o signed int array with query kmers (by position)
|   - lenQrySI:
|     o number kmers in qryKmerArySI
|   - refKmerArySI:
|     o signed int array with reference kmers to chain
|   - refIndexArySI:
|     o signed int array with refKmerArySI index's sorted
|       by kmer
|     o create with fun28 sortKmerIndex_mapRead
|   - lenRefSI:
|     o number kmers in refKmerArySI
|   - refStartSI:
|     o first refrence base in window
| Output:
|   - Modifies:
|     o refKmerArySI, qryKmerArySI, and chainLenArySI in
|       chainsSTPtr to have chain starting coordinates and
|       length
|     o scoreArySI and nextArySI in chainsSTPtr to have
|       all detected kmer positions filled with 0's
|     o indexArySI in chainsSTPtr to have index or -1 for
|       each query kmer start position in qryKmerArySI
|     o lenChainSI and lenIndexSI to have length of
|       indexArySI (num ref kmers) and number of chains
|     o sizeChainSI and sizeIndexSI in chainsSTPtr if
|       needed to rellocate memory for chain/index arrays
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|   - Note:
|     o sorting order of chains is query coordinate and
|       then reference coordinate
\-------------------------------------------------------*/
signed char
kmerChain_mapRead(
   struct chains_mapRead *chainsSTPtr,/*gets kmer chains*/
   signed int minChainLenSI,          /*min chain length*/
   unsigned int kmerUC,               /*length of 1 kmer*/
   signed int *qryKmerArySI,          /*query kmers*/
   signed int lenQrySI,               /*query length*/
   signed int *refKmerArySI,          /*reference kmers*/
   signed int *refIndexArySI,
     /*sorted refernce kmer index's (fun28)*/
   signed int lenRefSI,               /*number ref kmers*/
   signed int refStartSI              /*first ref base*/
);

/*-------------------------------------------------------\
| Fun33: mergeChains_mapRead
|   - merges kmer chains from fun31 (kmerChain) into
|     the chain the maximizes base count
|   - not optimal solution, but also reduce complexity
|     by not counting for gaps
| Input:
|   - chainsSTPtr:
|     o chains_mapRead struct pionter with chains to
|       merge and to find "best" chain in
|   - lenKmerUC:
|     o size of one kmer
|   - gapSI:
|     o gap penalty to apply for indels between chains
|   - matchSI:
|     o score for match
|   - maxLenSI:
|     o maximum length between chains before merging
| Output:
|   - Modifies:
|     o nextArySI in chains_mapRead to have best scoring
|       trail for each chain
|     o scoreArySI to have most bases possible in a chain
|   - Returns:
|     o index of best merged chain
|     o def_noChains_mapRead if no chains in array
\-------------------------------------------------------*/ 
signed int
mergeChains_mapRead(
   struct chains_mapRead *chainsSTPtr, /*chains to merge*/
   unsigned char lenKmerUC,            /*size of 1 kmer*/
   signed char gapSI,                  /*score for gaps*/
   signed char matchSI,                /*score for match*/
   signed int maxLenSI          /*length between chains*/
);

/*-------------------------------------------------------\
| Fun34: mergeToSam_mapRead
|   - convert output from mergeChains_map read to sam
|     file entry (mapq is set to 0)
| Input:
|   - chainsSTPtr:
|     o chains_mapRead struct pointer with chains to use
|       in making the alignment
|   - indexSI:
|     o starting index of merged chain to align
|   - alnEndsBl:
|     o align up to 1.5 chain lengths of the softmasked
|       ends of reads (0 = do not align)
|   - qrySTPtr:
|     o seqST struct with query sequence
|   - refSTPtr:
|     o seqST struct pointer with reference sequence
|   - lenKmerUC:
|     o length of one kmer in reference/query
|   - chainLenSI:
|     o minimum length of a chain
|     o used to align ends before chains, alignment size
|       is 1.5 chains (for chainLenSI = 20, is 30 bases)
|       * 1.5 is off the top of my head
|   - settings:
|     o settings for alignment
|   - matrixSTPtr:
|     o directional matrix to use in alignment step
|   - samSTPtr:
|     o samEntry struct pointer to add alignment to
| Output:
|   - Modifies:
|     o samSTPtr to have alignment; everything in samSTPtr
|       is modified
|     o matrixSTPtr (nothing usefull, here to reduce
|       memory allocations in alignment step)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
\-------------------------------------------------------*/
signed char
mergeToSam_mapRead(
   struct chains_mapRead *chainsSTPtr,/*chains to use*/
   signed int indexSI,            /*chain index to start*/
   signed char alnEndsBl,      /*1: map end masks*/
   struct seqST *qrySTPtr,     /*query sequence*/
   struct seqST *refSTPtr,     /*reference sequence*/
   unsigned char lenKmerUC,    /*length of one kmer*/
   signed int chainLenSI,      /*min length of one chain*/
   struct alnSet *settings,      /*aligment settings*/
   struct dirMatrix *matrixSTPtr,/*for alignment step*/
   struct samEntry *samSTPtr     /*gets final alignment*/
);

/*-------------------------------------------------------\
| Fun35: scoreSubAln_mapRead
|   - scores the aligment in samEntry and finds best
|     sub-aligment
| Input:
|   - samSTPtr:
|     o samEntry struct pointer with alignment to score
|     o do not trim
|   - subBl:
|     o modifiy alignment to best sub-alignment
|       * 1: look for best sub-alignment
|       * 0: do not look for sub-alignments
|   - qrySTPtr:
|     o seqST struct pointer with query sequence used in
|       alignment
|     o seqStr must be as lookup index (seqToIndex_alnSet)
|   - refSTPtr:
|     o seqST struct pointer with reference sequence used
|       in alignment
|     o seqStr must be as lookup index (seqToIndex_alnSet)
|   - alnSetSTPTr:
|     o alnSet struct pointer with scoring settings
|   - maxScoreSLPtr:
|     o signed long pointer to hold maximum score possible
|       for the aligned part of the query
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o extraStr in samSTPtr to have score "AS:i:score"
|       entry
|     o cigTypeStr and cigArySI in samSTPtr to be setup
|       for best sub-alignment (if subBl is 1)
|     o lenCigUI in samSTPtr to be new cigar length if
|       looking for best sub-alignment
|     o maxScoreSLPtr to have maximum possible score for
|       returned alignment (exlcudes unaligned bases)
|     o errSCPtr:
|       * 0 for no errors
|       * def_noScore_mapRead if no sub-score greater
|         than 0 (unlikley)
|       * def_memErr_mapRead for memory errors
|   - Returns:
|     o score for best sub-alignment or if subSl is 0, for
|       the entire samEntry
\-------------------------------------------------------*/
signed long
scoreSub_mapRead(
   struct samEntry *samSTPtr, /*alignment to score*/
   signed char subBl,         /*1:look for sub-alingment*/
   struct seqST *qrySTPtr,    /*query sequence*/
   struct seqST *refSTPtr,    /*reference sequence*/
   struct alnSet *alnSetSTPtr,/*scoring settings*/
   signed long *maxScoreSLPtr,/*gets max possible score*/
   signed char *errSCPtr      /*errors*/
);

/*-------------------------------------------------------\
| Fun36: chainToAln_mapRead
|   - converts chain in chains_mapRead struct to alignment
| Input:
|   - samSTPtr:
|     o samEntry struct pointer to add alignment to
|   - chainsSTPtr:
|     o chains_mapRead struct pointer with chains to use
|       in making the alignment
|   - indexSI:
|     o starting index of merged chain to align
|   - qrySTPtr:
|     o seqST struct with query sequence
|   - refSTPtr:
|     o seqST struct pointer with reference sequence
|   - setSTPtr:
|     o set_mapRead struct pointer with settings to use
|   - lenKmerUC:
|     o length of one kmer in reference/query
|   - chainLenSI:
|     o minimum length of a chain
|     o used to align ends before chains, alignment size
|       is 1.5 chains (for chainLenSI = 20, is 30 bases)
|       * 1.5 is off the top of my head
|   - matrixSTPtr:
|     o directional matrix to use in alignment step
|   - errSCPtr:
|     o signed char pointer to hold errors
| Output:
|   - Modifies:
|     o samSTPtr to have alignment (if could align)
|     o uses matrixSTPtr for aligment step, is modified,
|       and resized, but nothing of relavence
|     o sets errSCPtr
|       * 0 for no errors
|       * def_memErr_mapRead for memory errors
|       * def_noChains_mapRead if merged chains at index
|         did not have enough coverage
|       * def_noAln_mapRead if aligment was under
|         minimum percent query length or score
|   - Returns:
|     o 0 for no errors
|     o def_memErr_mapRead for memory errors
|     o def_noChains_mapRead if merged chains at index did
|       not have enough coverage
|     o def_noAln_mapRead if aligment was under
|       minimum percent query length or score
\-------------------------------------------------------*/
signed long
chainToAln_mapRead(
   struct samEntry *samSTPtr,    /*gets alignment*/
   struct chains_mapRead *chainsSTPtr, /*chains to align*/
   signed int indexSI,           /*chain index to align*/
   struct seqST *qrySTPtr,       /*query sequence*/
   struct seqST *refSTPtr,       /*reference sequence*/
   struct set_mapRead *setSTPtr, /*alignment settings*/
   unsigned char lenKmerUC,      /*lenght of 1 kemr*/
   signed int chainLenSI,        /*min length of a chain*/
   struct dirMatrix *matrixSTPtr,/*for alignment*/
   signed char *errSCPtr         /*gets errors*/
);

/*-------------------------------------------------------\
| Fun37: align_mapRead
|   - maps read to reference
| Input:
|   - qrySTPtr:
|     o seqST struct pointer with query sequence and id
|   - refSTPtr:
|     o ref_mapRead struct pointer with reference
|       sequence, settings, and kmer arrays for alignment
|     o setup with init_ref_mapRead, setup_ref_mapRead,
|       and addRef_ref_mapRead
|   - mapIndexSIPtr:
|     o signed int pointer to hold index of best chain
|       (used to make mapping)
|   - alnSTPtr:
|     o aln_mapRead struct pionter with allocated memory to
|       use/allocate/resize in mapping steps
|   - setSTPtr:
|     o set_mapRead structure with non-reference specific
|   - samSTPtr:
|     o pointer to samEntry struct to hold best alignment
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o samSTPtr to have highest scoring alignment
|     o seqStr in qrySTPtr to be lookup index's, you can
|       undo this with indexToSeq_alnSet(qrySTPtr->seqStr)
|     o mapIndexSLPtr to have index of chain used to find
|       the best alignment
|     o kmerArySI in alnSTPtr to have query kmer sequence
|     o lenSI in alnSTPtr to be length of kmerArySI
|     o sizeSI in alnSTPtr if kmerArySI is resized
|     o chainsSTPtr alnSTPtr to have found chains
|     o errSCPtr to have error message
|       * 0 for no errors
|       * def_memErr_mapRead for memory errors
|       * def_noChains_mapRead if chould not find any
|         chains
|       * def_noAln_mapRead if aligment was under
|         minimum percent query length or score
|     o matrixSTPtr in alnSTPtr to nothing usefull. Here
|       to avoid memory allocations in alignment step
|   - Returns:
|     o score for primary alignment
\-------------------------------------------------------*/
signed long
align_mapRead(
   struct seqST *qrySTPtr,      /*has query sequence*/
   struct ref_mapRead *refSTPtr,/*ref kmers/settings*/
   signed int *mapIndexSIPtr,   /*index of best chain*/
   struct samEntry *samSTPtr,   /*gets alignments*/
   struct aln_mapRead *alnSTPtr,/*memory for query*/
   struct set_mapRead *setSTPtr,/*alignment settings*/
   signed char *errSCPtr        /*gets errors*/
);

#endif

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconveint / not possible, this code is under the
:   MIT license
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy, modify, publish, use, compile,
:   sell, or distribute this software, either in source
:   code form or as a compiled binary, for any purpose,
:   commercial or non-commercial, and by any means.
: 
: In jurisdictions that recognize copyright laws, the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
:   CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information, please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2025 jeremyButtler
: 
: Permission is hereby granted, free of charge, to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"), to
:   deal in the Software without restriction, including
:   without limitation the rights to use, copy, modify,
:   merge, publish, distribute, sublicense, and/or sell
:   copies of the Software, and to permit persons to whom
:   the Software is furnished to do so, subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
:   AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
:   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
