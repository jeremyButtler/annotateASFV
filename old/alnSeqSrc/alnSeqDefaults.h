/*#########################################################
# Name: alnSeqDefaults.h
# Use:
#  o Holds the default settings and global definitions for
#    find alignSeq
# Libraries:
# C Standard Libraries:
#########################################################*/

#ifndef ALNSEQDEFAULTS_H
#define ALNSEQDEFAULTS_H

/*Gereral find co-infections settings*/
#define defVersion 20231025  /*Version number for alnSeq*/

// Aligment method being used
#define defUseNeedle 1  /*do Needleman-Wunsch alignment*/
#define defUseWater 0   /*Do a waterman smith alignment*/
#define defUseHirsch 0  /*Use a Hirschberg alignment*/
#define defUseMemWater 0 /*Use a memory efficent water*/
#define defJustScoresBl 0 /*memory efficent water only*/
   /*Only print out coordiantes and scores for memWater*/

/*Output file formating (only change defFormat)*/
#define defExpandCig 0 /*Expanded cigar format flag*/
#define defEMBOSS 1    /*EMBOSS format flag*/
#define defClustal 2   /*Clustal format flag*/
#define defFasta 3     /*Fasta format flag*/
#define defFormat defExpandCig
   /*Change defFormat to change the default output
   ` The current default format is expand cigar
   */

/*Doing printing out more than one aligment*/
#define defQueryRefScan 0
   /* 1: Recored the best score for each reference and
   `     query base in in a Waterman Smith alnigment
   */
#define defMinScore 1000
    /*Min score to keep an non-best alignment when doing
    ` an mutli-entry Waterman Smith alignment
    */


/*Output settings*/
#define defPAln 1     /*1: Print out the entire alignment*/
#define defPPos 0     /*1: Print base positions*/
#define defLineWrap 59 /*Max chars per line in out file*/
   /*Minimun is 40 characters*/

/*Alignment matrix movements*/
/*Do not change these values*/
#define defMvStop 0    /*Stop*/
#define defMvDel 1     /*Move left (deletion)*/
#define defMvIns 2     /*Move up (insertion)*/
#define defMvSnp 3     /*Move on a diagnol (snp/match)*/

/*Matrix filling settings*/

/* Settings for which direction to keep when their are
`  multiple equal directions.
   Order is AlwaysKeep->secondBest->lastOption
*/
#define defSnpInsDel 0
#define defSnpDelIns 1
#define defInsSnpDel 2
#define defInsDelSnp 3
#define defDelSnpIns 4
#define defDelInsSnp 5

#define defBestDir defInsDelSnp

/*Scoring variables*/
#define defGapOpen -10  /*Penalty for starting indel*/

#if !defined NOGAPOPEN
   #define defGapExtend -1 /*penalty for extending indel*/
#else
   /*Need a harsher score when I have no gap open penalty*/
   #define defGapExtend -6 /*penalty for extending indel*/
#endif

/* Scoring matrix is EDNAFULL or something close to it. I
`   think I have the anonymous bases with the correct
`  scores.
*/

/*Scoring matrix name*/
static char defMatrixNameStr[] = "EDNAFULL";

/*Alignment scoring matrix (non-anonymous)*/
#define defAToA 5     /*Score for query A to ref A*/
#define defAToT -4    /*Score for query A to ref T*/
#define defAToU -4    /*Score for query A to ref T*/
#define defAToG -4    /*Score for query A to ref G*/
#define defAToC -4    /*Score for query A to ref C*/

#define defTToA -4    /*Score for query T to ref A*/
#define defTToT 5     /*Score for query T to ref T*/
#define defTToU 5     /*Score for query T to ref T*/
#define defTToG -4    /*Score for query T to ref G*/
#define defTToC -4    /*Score for query T to ref C*/

#define defGToA -4    /*Score for query G to ref A*/
#define defGToT -4    /*Score for query G to ref T*/
#define defGToU -4    /*Score for query G to ref T*/
#define defGToG 5     /*Score for query G to ref G*/
#define defGToC -4    /*Score for query G to ref C*/

#define defCToA -4    /*Score for query C to ref A*/
#define defCToT -4    /*Score for query C to ref T*/
#define defCToU -4    /*Score for query C to ref T*/
#define defCToG -4    /*Score for query C to ref G*/
#define defCToC 5     /*Score for query C to ref C*/

/*Alignment scoring matrix anonymous (A)*/
#define defAToW 1     /*Score for query A to ref W (AT)*/
#define defAToS -4    /*Score for query A to ref S (CG)*/
#define defAToM 1     /*Score for query A to ref M (AC)*/
#define defAToK -4    /*Score for query A to ref K (GT)*/
#define defAToR 1     /*Score for query A to ref R (AG)*/
#define defAToY -4    /*Score for query A to ref Y (CT)*/
#define defAToB -4    /*Score for query A to ref B (CGT)*/
#define defAToD -1    /*Score for query A to ref D (AGT)*/
#define defAToH -1    /*Score for query A to ref H (ACT)*/
#define defAToV -1    /*Score for query A to ref V (ACG)*/
#define defAToN -2    /*Score for query A to ref N (ACGT)*/
#define defAToX -2    /*Score for query A to ref X (ACGT)*/

#define defWToA 1     /*Score for query W to ref A (AT)*/
#define defSToA -4    /*Score for query S to ref A (CG)*/
#define defMToA 1     /*Score for query M to ref A (AC)*/
#define defKToA -4    /*Score for query K to ref A (GT)*/
#define defRToA 1     /*Score for query R to ref A (AG)*/
#define defYToA -4    /*Score for query Y to ref A (CT)*/
#define defBToA -4    /*Score for query B to ref A (CGT)*/
#define defDToA -1    /*Score for query D to ref A (AGT)*/
#define defHToA -1    /*Score for query H to ref A (ACT)*/
#define defVToA -1    /*Score for query V to ref A (ACG)*/
#define defNToA -2    /*Score for query N to ref A (ACGT)*/
#define defXToA -2    /*Score for query X to ref A (ACGT)*/

/*Alignment scoring matrix anonymous (C)*/
#define defCToW -4    /*Score for query C to ref W (AT)*/
#define defCToS 1     /*Score for query C to ref S (CG)*/
#define defCToM 1     /*Score for query C to ref M (AC)*/
#define defCToK -4    /*Score for query C to ref K (GT)*/
#define defCToR -4    /*Score for query C to ref R (AG)*/
#define defCToY 1     /*Score for query C to ref Y (CT)*/
#define defCToB -1    /*Score for query C to ref B (CGT)*/
#define defCToD -4    /*Score for query C to ref D (AGT)*/
#define defCToH -1    /*Score for query C to ref H (ACT)*/
#define defCToV -1    /*Score for query C to ref V (ACG)*/
#define defCToN -2    /*Score for query C to ref N (ACGT)*/
#define defCToX -2    /*Score for query C to ref X (ACGT)*/

#define defWToC -4    /*Score for query W to ref C (AT)*/
#define defSToC 1     /*Score for query S to ref C (CG)*/
#define defMToC 1     /*Score for query M to ref C (AC)*/
#define defKToC -4    /*Score for query K to ref C (GT)*/
#define defRToC -4    /*Score for query R to ref C (AG)*/
#define defYToC 1     /*Score for query Y to ref C (CT)*/
#define defBToC -1    /*Score for query B to ref C (CGT)*/
#define defDToC -4    /*Score for query D to ref C (AGT)*/
#define defHToC -1    /*Score for query H to ref C (ACT)*/
#define defVToC -1    /*Score for query V to ref C (ACG)*/
#define defNToC -2    /*Score for query N to ref C (ACGT)*/
#define defXToC -2    /*Score for query X to ref C (ACGT)*/

/*Alignment scoring matrix anonymous (G)*/
#define defGToW -4    /*Score for query G to ref W (AT)*/
#define defGToS 1     /*Score for query G to ref S (CG)*/
#define defGToM -4    /*Score for query G to ref M (AC)*/
#define defGToK 1     /*Score for query G to ref K (GT)*/
#define defGToR 1     /*Score for query G to ref R (AG)*/
#define defGToY -4    /*Score for query G to ref Y (CT)*/
#define defGToB -1    /*Score for query G to ref B (CGT)*/
#define defGToD -1    /*Score for query G to ref D (AGT)*/
#define defGToH -4    /*Score for query G to ref H (ACT)*/
#define defGToV -1    /*Score for query G to ref V (ACG)*/
#define defGToN -2    /*Score for query G to ref N (ACGT)*/
#define defGToX -2    /*Score for query G to ref X (ACGT)*/

#define defWToG -4    /*Score for query G to ref W (AT)*/
#define defSToG 1     /*Score for query G to ref S (CG)*/
#define defMToG -4    /*Score for query G to ref M (AC)*/
#define defKToG 1     /*Score for query G to ref K (GT)*/
#define defRToG 1     /*Score for query G to ref R (AG)*/
#define defYToG -4    /*Score for query G to ref Y (CT)*/
#define defBToG -1    /*Score for query G to ref B (CGT)*/
#define defDToG -1    /*Score for query G to ref D (AGT)*/
#define defHToG -4    /*Score for query G to ref H (ACT)*/
#define defVToG -1    /*Score for query G to ref V (ACG)*/
#define defNToG -2    /*Score for query G to ref N (ACGT)*/
#define defXToG -2    /*Score for query G to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (T)*/
#define defTToW 1     /*Score for query T to ref W (AT)*/
#define defTToS -4    /*Score for query T to ref S (CG)*/
#define defTToM -4    /*Score for query T to ref M (AC)*/
#define defTToK 1     /*Score for query T to ref K (GT)*/
#define defTToR -4    /*Score for query T to ref R (AG)*/
#define defTToY 1     /*Score for query T to ref Y (CT)*/
#define defTToB -1    /*Score for query T to ref B (CGT)*/
#define defTToD -1    /*Score for query T to ref D (AGT)*/
#define defTToH -1    /*Score for query T to ref H (ACT)*/
#define defTToV -4    /*Score for query T to ref V (ACG)*/
#define defTToN -2    /*Score for query T to ref N (ACGT)*/
#define defTToX -2    /*Score for query T to ref X (ACGT)*/

#define defWToT 1     /*Score for query T to ref W (AT)*/
#define defSToT -4    /*Score for query T to ref S (CG)*/
#define defMToT -4    /*Score for query T to ref M (AC)*/
#define defKToT 1     /*Score for query T to ref K (GT)*/
#define defRToT -4    /*Score for query T to ref R (AG)*/
#define defYToT -4    /*Score for query T to ref Y (CT)*/
#define defBToT -1    /*Score for query T to ref B (CGT)*/
#define defDToT -1    /*Score for query T to ref D (AGT)*/
#define defHToT -1    /*Score for query T to ref H (ACT)*/
#define defVToT -4    /*Score for query T to ref V (ACG)*/
#define defNToT -2    /*Score for query T to ref N (ACGT)*/
#define defXToT -2    /*Score for query T to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (U)*/
#define defUToW 1     /*Score for query T to ref W (AT)*/
#define defUToS -4    /*Score for query T to ref S (CG)*/
#define defUToM -4    /*Score for query T to ref M (AC)*/
#define defUToK 1     /*Score for query T to ref K (GT)*/
#define defUToR -4    /*Score for query T to ref R (AG)*/
#define defUToY 1     /*Score for query T to ref Y (CT)*/
#define defUToB -1    /*Score for query T to ref B (CGT)*/
#define defUToD -1    /*Score for query T to ref D (AGT)*/
#define defUToH -1    /*Score for query T to ref H (ACT)*/
#define defUToV -4    /*Score for query T to ref V (ACG)*/
#define defUToN -2    /*Score for query T to ref N (ACGT)*/
#define defUToX -2    /*Score for query T to ref X (ACGT)*/

#define defWToU 1     /*Score for query T to ref W (AT)*/
#define defSToU -4    /*Score for query T to ref S (CG)*/
#define defMToU -4    /*Score for query T to ref M (AC)*/
#define defKToU 1     /*Score for query T to ref K (GT)*/
#define defRToU -4    /*Score for query T to ref R (AG)*/
#define defYToU -4    /*Score for query T to ref Y (CT)*/
#define defBToU -1    /*Score for query T to ref B (CGT)*/
#define defDToU -1    /*Score for query T to ref D (AGT)*/
#define defHToU -1    /*Score for query T to ref H (ACT)*/
#define defVToU -4    /*Score for query T to ref V (ACG)*/
#define defNToU -2    /*Score for query T to ref N (ACGT)*/
#define defXToU -2    /*Score for query T to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (W)*/
#define defWToW -1    /*Score for query W to ref W (AT)*/
#define defWToS -4    /*Score for query W to ref S (CG)*/
#define defWToM -1    /*Score for query W to ref M (AC)*/
#define defWToK -1    /*Score for query W to ref K (GT)*/
#define defWToR -1    /*Score for query W to ref R (AG)*/
#define defWToY -1    /*Score for query W to ref Y (CT)*/
#define defWToB -1    /*Score for query W to ref B (CGT)*/
#define defWToD -1    /*Score for query W to ref D (AGT)*/
#define defWToH -1    /*Score for query W to ref H (ACT)*/
#define defWToV -1    /*Score for query W to ref V (ACG)*/
#define defWToN -1    /*Score for query W to ref N (ACGT)*/
#define defWToX -1    /*Score for query W to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (S)*/
#define defSToW -4    /*Score for query S to ref W (AT)*/
#define defSToS -1    /*Score for query S to ref S (CG)*/
#define defSToM -2    /*Score for query S to ref M (AC)*/
#define defSToK -2    /*Score for query S to ref K (GT)*/
#define defSToR -2    /*Score for query S to ref R (AG)*/
#define defSToY -2    /*Score for query S to ref Y (CT)*/
#define defSToB -1    /*Score for query S to ref B (CGT)*/
#define defSToD -3    /*Score for query S to ref D (AGT)*/
#define defSToH -3    /*Score for query S to ref H (ACT)*/
#define defSToV -1    /*Score for query S to ref V (ACG)*/
#define defSToN -1    /*Score for query S to ref N (ACGT)*/
#define defSToX -1    /*Score for query S to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (M)*/
#define defMToW -2    /*Score for query M to ref W (AT)*/
#define defMToS -2    /*Score for query M to ref S (CG)*/
#define defMToM -1    /*Score for query M to ref M (AC)*/
#define defMToK -4    /*Score for query M to ref K (GT)*/
#define defMToR -2    /*Score for query M to ref R (AG)*/
#define defMToY -2    /*Score for query M to ref Y (CT)*/
#define defMToB -3    /*Score for query M to ref B (CGT)*/
#define defMToD -3    /*Score for query M to ref D (AGT)*/
#define defMToH -1    /*Score for query M to ref H (ACT)*/
#define defMToV -1    /*Score for query M to ref V (ACG)*/
#define defMToN -1    /*Score for query M to ref N (ACGT)*/
#define defMToX -1    /*Score for query M to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (K)*/
#define defKToW -2    /*Score for query K to ref W (AT)*/
#define defKToS -2    /*Score for query K to ref S (CG)*/
#define defKToM -4    /*Score for query K to ref M (AC)*/
#define defKToK -1    /*Score for query K to ref K (GT)*/
#define defKToR -2    /*Score for query K to ref R (AG)*/
#define defKToY -2    /*Score for query K to ref Y (CT)*/
#define defKToB -1    /*Score for query K to ref B (CGT)*/
#define defKToD -1    /*Score for query K to ref D (AGT)*/
#define defKToH -3    /*Score for query K to ref H (ACT)*/
#define defKToV -3    /*Score for query K to ref V (ACG)*/
#define defKToN -1    /*Score for query K to ref N (ACGT)*/
#define defKToX -1    /*Score for query K to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (R)*/
#define defRToW -1    /*Score for query R to ref W (AT)*/
#define defRToS -1    /*Score for query R to ref S (CG)*/
#define defRToM -1    /*Score for query R to ref M (AC)*/
#define defRToK -2    /*Score for query R to ref K (GT)*/
#define defRToR -1    /*Score for query R to ref R (AG)*/
#define defRToY -4    /*Score for query R to ref Y (CT)*/
#define defRToB -3    /*Score for query R to ref B (CGT)*/
#define defRToD -1    /*Score for query R to ref D (AGT)*/
#define defRToH -1    /*Score for query R to ref H (ACT)*/
#define defRToV -3    /*Score for query R to ref V (ACG)*/
#define defRToN -1    /*Score for query R to ref N (ACGT)*/
#define defRToX -1    /*Score for query R to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (Y)*/
#define defYToW -2    /*Penalty for query Y to ref W (AT)*/
#define defYToS -2    /*Score for query Y to ref S (CG)*/
#define defYToM -2    /*Score for query Y to ref M (AC)*/
#define defYToK -2    /*Score for query Y to ref K (GT)*/
#define defYToR -4    /*Score for query Y to ref R (AG)*/
#define defYToY -1    /*Score for query Y to ref Y (CT)*/
#define defYToB -2    /*Score for query Y to ref B (CGT)*/
#define defYToD -3    /*Score for query Y to ref D (AGT)*/
#define defYToH -1    /*Score for query Y to ref H (ACT)*/
#define defYToV -3    /*Score for query Y to ref V (ACG)*/
#define defYToN -1    /*Score for query Y to ref N (ACGT)*/
#define defYToX -1    /*Score for query Y to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (B)*/
#define defBToW -3    /*Score for query B to ref W (AT)*/
#define defBToS -1    /*Score for query B to ref S (CG)*/
#define defBToM -3    /*Score for query B to ref M (AC)*/
#define defBToK -1    /*Score for query B to ref K (GT)*/
#define defBToR -3    /*Score for query B to ref R (AG)*/
#define defBToY -1    /*Score for query B to ref Y (CT)*/
#define defBToB -1    /*Score for query B to ref B (CGT)*/
#define defBToD -2    /*Score for query B to ref D (AGT)*/
#define defBToH -2    /*Score for query B to ref H (ACT)*/
#define defBToV -2    /*Score for query B to ref V (ACG)*/
#define defBToN -1    /*Score for query B to ref N (ACGT)*/
#define defBToX -1    /*Score for query B to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (D)*/
#define defDToW -1    /*Score for query D to ref W (AT)*/
#define defDToS -3    /*Score for query D to ref S (CG)*/
#define defDToM -3    /*Score for query D to ref M (AC)*/
#define defDToK -1    /*Score for query D to ref K (GT)*/
#define defDToR -1    /*Score for query D to ref R (AG)*/
#define defDToY -3    /*Score for query D to ref Y (CT)*/
#define defDToB -2    /*Score for query D to ref B (CGT)*/
#define defDToD -1    /*Score for query D to ref D (AGT)*/
#define defDToH -2    /*Score for query D to ref H (ACT)*/
#define defDToV -2    /*Score for query D to ref V (ACG)*/
#define defDToN -1    /*Score for query D to ref N (ACGT)*/
#define defDToX -1    /*Score for query D to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (H)*/
#define defHToW -1    /*Score for query H to ref W (AT)*/
#define defHToS -3    /*Score for query H to ref S (CG)*/
#define defHToM -1    /*Score for query H to ref M (AC)*/
#define defHToK -3    /*Score for query H to ref K (GT)*/
#define defHToR -3    /*Score for query H to ref R (AG)*/
#define defHToY -1    /*Score for query H to ref Y (CT)*/
#define defHToB -2    /*Score for query H to ref B (CGT)*/
#define defHToD -2    /*Score for query H to ref D (AGT)*/
#define defHToH -1    /*Score for query H to ref H (ACT)*/
#define defHToV -2    /*Score for query H to ref V (ACG)*/
#define defHToN -1    /*Score for query H to ref N (ACGT)*/
#define defHToX -1    /*Score for query H to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (V)*/
#define defVToW -3    /*Score for query V to ref W (AT)*/
#define defVToS -1    /*Score for query V to ref S (CG)*/
#define defVToM -1    /*Score for query V to ref M (AC)*/
#define defVToK -3    /*Score for query V to ref K (GT)*/
#define defVToR -1    /*Score for query V to ref R (AG)*/
#define defVToY -3    /*Score for query V to ref Y (CT)*/
#define defVToB -2    /*Score for query V to ref B (CGT)*/
#define defVToD -2    /*Score for query V to ref D (AGT)*/
#define defVToH -2    /*Score for query V to ref H (ACT)*/
#define defVToV -1    /*Score for query V to ref V (ACG)*/
#define defVToN -1    /*Score for query V to ref N (ACGT)*/
#define defVToX -1    /*Score for query V to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (N)*/
#define defNToW -2    /*Score for query N to ref W (AT)*/
#define defNToS -2    /*Score for query N to ref S (CG)*/
#define defNToM -2    /*Score for query N to ref M (AC)*/
#define defNToK -2    /*Score for query N to ref K (GT)*/
#define defNToR -2    /*Score for query N to ref R (AG)*/
#define defNToY -2    /*Score for query N to ref Y (CT)*/
#define defNToB -2    /*Score for query N to ref B (CGT)*/
#define defNToD -2    /*Score for query N to ref D (AGT)*/
#define defNToH -2    /*Score for query N to ref H (ACT)*/
#define defNToV -2    /*Score for query N to ref V (ACG)*/
#define defNToN -1    /*Score for query N to ref N (ACGT)*/
#define defNToX -1    /*Score for query N to ref X (ACGT)*/

/*Alignment scoring matrix anonymous (X)*/
#define defXToW -1    /*Score for query X to ref W (AT)*/
#define defXToS -1    /*Score for query X to ref S (CG)*/
#define defXToM -1    /*Score for query X to ref M (AC)*/
#define defXToK -1    /*Score for query X to ref K (GT)*/
#define defXToR -1    /*Score for query X to ref R (AG)*/
#define defXToY -1    /*Score for query X to ref Y (CT)*/
#define defXToB -1    /*Score for query X to ref B (CGT)*/
#define defXToD -1    /*Score for query X to ref D (AGT)*/
#define defXToH -1    /*Score for query X to ref H (ACT)*/
#define defXToV -1    /*Score for query X to ref V (ACG)*/
#define defXToN -1    /*Score for query X to ref N (ACGT)*/
#define defXToX -1    /*Score for query X to ref X (ACGT)*/

#endif
