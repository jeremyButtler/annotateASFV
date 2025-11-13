/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' seqStats SOF: Start Of File
'   - get stats for fx (fastq/fasta) files
'   o header:
'     - guards and forward declarations
'   o .c fun01: ln_seqStats
'     - natureal log value without using math.h (requies
'       -lm flag)
'     - TODO: fix issues with low values
'   o .c fun02: log10_seqStats
'     - log10 value without using math.h (requies -lm
'       flag)
'   o .c fun03: pow10_seqStats
'     - rases 10 to the input exponant
'   o .c fun04: qhistToMedQ_seqStats
'     - converts a q-score histogram to medain q-score
'   o .c fun05: qhistToMeanQ_seqStats
'     - converts a q-score histogram to mean q-score
'   o fun06: meanMedQ_seqStats
'     - gets the mean and median q-score for a sequence
'   o fun07: ontMeanMedQ_seqStats
'     - gets the mean and median q-score for a sequence
'       using ONT's method
'   o fun08: pReadStats_seqStats
'     - prints stats for a single read
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards and forward declarations
\-------------------------------------------------------*/

#ifndef SEQUENCE_STATS_H
#define SEQUENCE_STATS_H

struct seqST;

/*-------------------------------------------------------\
| Fun06: meanMedQ_seqStats
|   - gets the mean and median q-score for a sequence
| Input:
|   - qStr:
|     o c-string with q-score entry to get q-scores for
|   - seqLenSI:
|     o number of bases in the sequence (index 1)
|   - medianQFPtr:
|     o float pointer to get the median q-score
| Output:
|   - Modifies:
|     o medainQFPtr to have the median q-score
|   - Returns:
|     o mean q-score
\-------------------------------------------------------*/
float
meanMedQ_seqStats(
   signed char *qStr,
   signed int seqLenSI,
   float *medianQFPtr
);

/*-------------------------------------------------------\
| Fun07: ontMeanMedQ_seqStats
|   - gets the mean and median q-score for a sequence using
|     ONT's method
| Input:
|   - qStr:
|     o c-string with q-score entry to get q-scores for
|   - seqLenSI:
|     o number of bases in the sequence (index 1)
|   - medianQFPtr:
|     o float pointer to get the median q-score
|   - keepFirst60Bl:
|     o 1: keep the first 60 bases in the read
|     o 0: follow ONT and ingore first 60 bases
| Output:
|   - Modifies:
|     o medainQFPtr to have the median q-score
|   - Returns:
|     o mean q-score
\-------------------------------------------------------*/
float
ontMeanMedQ_seqStats(
   signed char *qStr,
   signed int seqLenSI,
   float *medianQFPtr,
   signed char keepFirst60Bl
);

/*-------------------------------------------------------\
| Fun08: pReadStats_seqStats
|   - prints stats for a single read
| Input:
|   - qStr:
|     o c-string with q-score entry to get q-scores for
|   - seqLenSI:
|     o number of bases in the sequence (index 1)
|   - idStr:
|     o c-string with id to print out
|   - pHeadBlPtr:
|     o 1: print header and set to 0
|     o 0: do not print header
|   - keepFirst60Bl:
|     o 1: for ONT mean keep the firt 60 bases
|     o 0: for ONT mean discard the first 60 bases
|   - statMedthodSC:
|     o 1: use ONT stat only
|     o 2: use mean/median (no error step) only
|     o 3: print both ONT and mean/median
|   - outFILE:
|     o FILE pionter to print stats to
| Output:
|   - Prints;
|     o read stats (and header if requested) to outFILE
|   - Modifies:
|     o pHeadBlPtr to be 0 (if is not 0)
\-------------------------------------------------------*/
void
pReadStats_seqStats(
   signed char *qStr,        /*quality score entry*/
   signed int seqLenSI,      /*bases in quality score*/
   signed char *idStr,       /*read id/name to print*/
   signed char *pHeadBlPtr,  /*1: print header*/
   signed char keepFirst60Bl,/*1: ONT keep first 60 nt*/
   signed char statMethodSC, /*1:ONT, 2:mean/med; 3:both*/
   void *outFILE             /*file to print to*/
);

#endif
