#include "hirschVect16.h"

/*--------------------------------------------------------\
| Name: HirschberVectI16
| Use:
|  - Runs a Hirschberg alignment with 16 bit vectors.
| Input:
|  - refST:
|    o seqStruct with the reference sequence and starting
|      and ending cooridinates
|  - qryST:
|    o qryStruct with the query sequence and starting and
|      ending cooridinates
|  - alnSet:
|    o Settings to run the alignment with.
| Output:
|  - Returns:
|    o A alignment structure with the alignment.
|    o 0 For memory errors
\--------------------------------------------------------*/
struct alnStruct * hirschVectI16(
  struct seqStruct *refST, /*Reference sequence to align*/
  struct seqStruct *qryST, /*Qeury sequence to align*/
    /* For refST and qryST, use seqStruct->offsetUL to set
    `  the starting point for the alignmnet and
    `  seqStruct->endAlnUL to set the ending point
    */
  struct alnSet *settings /*Settings to use for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-01 TOC: Hirschberg
   '  - Sets up for and calls the recursvie function to
   '    run a Hirschberg alignment
   '  o fun-01 sec-01:
   '    - Variable declerations
   '  o fun-01 sec-02:
   '    - Memory allocation (set up for Hirschberg)
   '  o fun-01 sec-03:
   '    - Run the hirschberg alignment
   '  o fun-01 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   int16_t lenRefS =
     (int16_t) (refST->endAlnUL - refST->offsetUL + 1);
     /*+1 to convert to index 1 (values are index 0)*/
   int16_t lenQryS =
     (int16_t) (qryST->endAlnUL - qryST->offsetUL + 1);
     /*+ 1 to convert to index 1 (values are index 0)*/

   int16_t *forwardScoreRowL = 0;
   int16_t *reverseScoreRowL = 0;
   int16_t *gapRowS = 0;

   struct alnStruct *alnST = 0;

   char *refAlnS = 0;
   char *qryAlnS = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-02:
   ^   - Memory allocation (set up for Hirschberg)
   ^   o fun-01 sec-02 sub-01:
   ^     - Initalize the ouput alignment structure 
   ^   o fun-01 sec-02 sub-02:
   ^     - Initalize the scoring rows
   ^   o fun-01 sec-02 sub-03:
   ^     - Initalize the direction rows
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-01:
   *  - Initalize the ouput alignment structure 
   \******************************************************/

   refAlnS = calloc(lenRefS + 1, sizeof(char));

   if(refAlnS == 0) return 0;

   #ifndef NOGAPOPEN
     gapRowS =
        calloc(lenRefS + defNum16BitElms,sizeof(int16_t));

     if(gapRowS == 0)
     { /*If I could not make another direction row*/
        free(refAlnS);
        refAlnS = 0;
        return 0;
     } /*If I could not make another direction row*/
   #endif

   qryAlnS = calloc(lenQryS + 1, sizeof(char));

   if(qryAlnS == 0)
   { /*If I could not make another direction row*/
      free(refAlnS);
      refAlnS = 0;

      #ifndef NOGAPOPEN
         free(gapRowS);
         gapRowS = 0;
      #endif

      return 0;
   } /*If I could not make another direction row*/

   /******************************************************\
   * Fun-01 Sec-02 Sub-02:
   *  - Initalize the scoring rows
   \******************************************************/

   /* I am using full length arrays to make the later
   `  steps eaiser. This takes more memory, but makes life
   `  nicer
   */

   forwardScoreRowL = 
      malloc(sizeof(int16_t) * (lenRefS+defNum16BitElms));

   if(forwardScoreRowL == 0)
   { /*If had a memory allocatoin error*/
     free(refAlnS);
     refAlnS = 0;
     free(qryAlnS);
     qryAlnS = 0;
 
     #ifndef NOGAPOPEN
        free(gapRowS);
        gapRowS = 0;
     #endif

     return 0;
   } /*If had a memory allocatoin error*/

   reverseScoreRowL =
      malloc(sizeof(int16_t) * (lenRefS+defNum16BitElms));

   if(reverseScoreRowL == 0)
   { /*If had a memory allocatoin error*/
     free(refAlnS);
     refAlnS = 0;
     free(qryAlnS);
     qryAlnS = 0;
 
     #ifndef NOGAPOPEN
        free(gapRowS);
        gapRowS = 0;
     #endif

     free(forwardScoreRowL);
     return 0;
   } /*If had a memory allocatoin error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-03:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Sening in offset values, because alignment array is
   ` sized to the alignmnet region
   */
   hirschVectI16Fun(
     refST->seqStr + refST->offsetUL,
     0,                /*1st reference base to align*/
     lenRefS,          /*Length of ref region to align*/
     qryST->seqStr + qryST->offsetUL,
     0,                /*1st query base to align*/
     lenQryS,          /*length of query target region*/
     forwardScoreRowL, /*For scoring*/
     reverseScoreRowL, /*For scoring*/
     refAlnS,          /*Holds the reference alignment*/
     qryAlnS,          /*Holds the query alignment*/
     gapRowS,          /*tells if using gap open/extend*/
     settings          /*Settings for the alignment*/
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-01 Sec-04:
   ^    - Clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   free(forwardScoreRowL);
   free(reverseScoreRowL);

   #ifndef  NOGAPOPEN
      free(gapRowS);
      gapRowS = 0;
   #endif

   alnST = hirschVectToAlnST(refST,qryST,refAlnS,qryAlnS);

   free(refAlnS);
   refAlnS = 0;
   free(qryAlnS);
   qryAlnS = 0;
 
   return alnST; /*Is 0 if twoBitAlnToAlnSt failed*/
} /*Hirschberg*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the output alignment
\--------------------------------------------------------*/
void hirschVectI16Fun(
  char *refSeqStr,  /*Reference sequence*/
  int16_t refStartS, /*index 0: 1st bast to algin in ref*/
  int16_t refLenS,   /*index 1 length of region to Align*/

  char *qrySeqStr,   /*Query sequence*/
  int16_t qryStartS,  /*index 0 Starting query base*/
  int16_t qryLenS,    /*index 1 Length of query region*/

  int16_t *forwardScoreRowS, /*Scores for forward row*/
  int16_t *reverseScoreRowS, /*Scores for reverse row*/

  char *refAlnST,     /*Holds output reference alignment*/
  char *qryAlnST,     /*Holds the output query alignment*/

  int16_t *gapRowS,   /*Tells if using gap open or extend*/
  struct alnSet *settings /*Settings to use for alignment*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-02 TOC: HirschbergFun
   '  - Does the recursive part of a Hirschberg alignment
   '  o fun-02 sec-01:
   '    - Variable declerations
   '  o fun-02 sec-02:
   '    - Check if on a leaf (final part of alignment
   '  o fun-02 sec-03:
   '    - Get scores
   '  o fun-02 sec-04:
   '    - Find the midpoint
   '  o fun-02 sec-05:
   '    - Run the next hirschberg alignment
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-01:
   ^    - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   uint8_t bitC = 0;
   int16_t forwardIndelColS = 0;
   int16_t reverseIndelColS = 0;
   int16_t midpointS = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-02:
   ^  - Check if on a leaf (final part of alignment
   ^  o fun-02 sec-02 sub-01:
   ^    - Handle cases were I have just insertions
   ^  o fun-02 sec-02 sub-02:
   ^    - Handle cases were I have just deletions
   ^  o fun-02 sec-02 sub-03:
   ^    - Handle cases were I have to align last ref base
   ^  o fun-02 sec-02 sub-04:
   ^  - Handle cases were I have to align last query base
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-01:
   *  - Handle cases were I have just insertions
   \******************************************************/

   if(refLenS == 0)
   { /*If all remaing bases are query insertions*/
     qryAlnST += qryStartS;
     vectAddGaps(qryAlnST, qryLenS, defGapFlag);
     return; /*Nothing else to do*/
   } /*If all remaing bases are query insertions*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-02:
   *  - Handle cases were I have just deletions
   \******************************************************/

   if(qryLenS == 0)
   { /*If all remaing bases are query deletions*/
     refAlnST += refStartS;
     vectAddGaps(refAlnST, refLenS, defGapFlag);
     return; /*Nothing else to do*/
   } /*If all remaing bases are query deletions*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-03:
   *  - Handle cases were I have to align last ref base
   \******************************************************/

   if(refLenS == 1)
   { /*If I have to align the last reference base*/

     if(qryLenS == 0)
     { /*If bases are aligned (one reference & one query)*/
        refAlnST += refStartS;
        *refAlnST = defGapFlag;
        return; /*Finished*/
     } /*If bases are aligned (one reference & one query)*/

     if(qryLenS == 1)
     { /*If bases are aligned (one reference & one query)*/
       qrySeqStr += qryStartS;
       refSeqStr += refStartS;

       if(checkIfBasesMatch(qrySeqStr, refSeqStr) == 1)
         bitC = defMatchFlag;

        else bitC = defSnpFlag;

        qryAlnST += qryStartS;
        *qryAlnST = bitC;

        refAlnST += refStartS;
        *refAlnST = bitC;

        return; /*Finished*/
     } /*If bases are aligned (one reference & one query)*/

     vectI16PosOneBase(
       *(refSeqStr + refStartS),/*ref base*/
       refStartS,              /*Position of ref base*/
       qrySeqStr,              /*first base of query*/
       qryStartS,              /*positoin of query*/
       qryLenS,                /*Length of the query*/
       refAlnST,                /*Array to hold alignment*/
       qryAlnST,                /*Array to hold alignment*/
       settings                 /*Has Scoring variables*/
     );

     return; /*This base is now aligned*/
   } /*If I have to align the last reference base*/

   /******************************************************\
   * Fun-02 Sec-02 Sub-04:
   *  - Handle cases were I have to align last query base
   \******************************************************/

   if(qryLenS == 1)
   { /*If I have to align the last query base*/
     if(refLenS == 0)
     { /*If bases are aligned (one reference & one query)*/
        qryAlnST += qryStartS;
        *qryAlnST = defGapFlag;
        return; /*Finished*/
     } /*If bases are aligned (one reference & one query)*/

     if(refLenS == 1)
     { /*If bases are aligned (one reference & one query)*/
       qrySeqStr += qryStartS;
       refSeqStr += refStartS;

       if(checkIfBasesMatch(qrySeqStr, refSeqStr) == 1)
         bitC = defMatchFlag;

        else bitC = defSnpFlag;

        qryAlnST += qryStartS;
        *qryAlnST = bitC;
        refAlnST += refStartS;
        *refAlnST = bitC;

        return; /*Finished*/
     } /*If bases are aligned (one reference & one query)*/

     vectI16PosOneBase(
       *(qrySeqStr + qryStartS),/*ref base*/
       qryStartS,              /*Position of ref base*/
       refSeqStr,              /*first base of reference*/
       refStartS,              /*positoin of query*/
       refLenS,                /*Length of the query*/
       qryAlnST,                /*Array to hold alignment*/
       refAlnST,                /*Array to hold alignment*/
       settings                 /*Has Scoring variables*/
     );

     return; /*Finshed aligning this query base*/
   } /*If I have to align the last query base*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-02 Sec-03:
   ^  - Get scores
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    forwardIndelColS = 
      scoreForwardHirsch(
        refSeqStr,       /*Entire reference sequence*/
        refStartS,       /*Starting base of ref target*/
        refLenS,         /*length of ref target region*/
        qrySeqStr,       /*Query seq with coordinates*/
        qryStartS,       /*Starting base of query target*/
        qryLenS / 2,     /*Length of query target region*/
        forwardScoreRowS, /*Array of scores to fill*/
        refAlnST,         /*direction row for gap extend*/
        settings          /*setttings to use*/
    ); /*Get the scores for the forward direction*/
    /*For -DNOGAPOPEN, refAlnST is ignored*/

    reverseIndelColS = 
      scoreReverseHirsch(
        refSeqStr,         /*Entire reference sequence*/
        refStartS,         /*Starting base of ref target*/
        refLenS,           /*length of ref target region*/
        qrySeqStr,         /*Query seq with coordinates*/
        qryStartS + (qryLenS / 2),/*new query start*/
        qryLenS - (qryLenS / 2),  /*New query length*/
        reverseScoreRowS,  /*Array of scores to fill*/
        gapRowS,            /*direction row for gap extend*/
        settings           /* setttings to use*/
      ); /*Get the scores for the reverse direction*/
      /*For -DNOGAPOPEN, gapRowS is ignored*/
      /* I can get away with queryLen/2 here, because 
      `  queryLen is index 1 and the function takes in
      `  an lenth 1 argument
      `  I made this section thread safe by using refAlnST
      `    for the forward score and gapRowS for the
      `    reverse row.
      */

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-04:
   ^   - Find the midpoint
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *(forwardScoreRowS + refStartS + refLenS - 1) +=
     reverseIndelColS;

   midpointS = refLenS;

   for(
     int16_t baseS = 0;
     baseS < refLenS - 1;
     ++baseS
   ) { /*Loop; add up all scores*/
     *(forwardScoreRowS + refStartS + baseS) += 
       *(reverseScoreRowS + refStartS + baseS + 1);
       /*The reverse row is already reversed*/

     if(
       *(forwardScoreRowS + refStartS + baseS) >
       *(forwardScoreRowS + refStartS + midpointS -1)
     ) midpointS = baseS + 1;
   } /*Loop; add up all scores*/

   *(reverseScoreRowS + refStartS) += forwardIndelColS;

   if(
     *(reverseScoreRowS + refStartS) >
     *(forwardScoreRowS + refStartS + midpointS - 1)
   ) midpointS = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^  Fun-02 Sec-05:
   ^    - Run the Hirschberg alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   hirschVectI16Fun(
     refSeqStr,      /*Full reference sequence*/
     refStartS,      /*Full reference sequence*/
     midpointS,      /*Length of new reference sequence*/
     qrySeqStr,      /*Full query sequence*/
     qryStartS,      /*Start of queyr target region*/
     qryLenS / 2,    /*Length of query target region*/
     forwardScoreRowS,/*For scoring*/
     reverseScoreRowS,/*Has last line of scores*/
     refAlnST,        /*direction row for gap extend*/
     qryAlnST,        /*Holds the alignment codes*/
     gapRowS,          /*For threadsafe scoreing*/
     settings         /*Settings for the alignment*/
   );

   hirschVectI16Fun(
     refSeqStr,              /*Full reference sequence*/
     refStartS + midpointS, /*New reference start*/
     refLenS - midpointS,   /*New reference end*/
     qrySeqStr,              /*Full query sequence*/
     qryStartS + (qryLenS / 2), /*New query start*/
     qryLenS - (qryLenS / 2),/*New query length*/
     forwardScoreRowS,        /*For scoring*/
     reverseScoreRowS,        /*Has last line of scores*/
     refAlnST,                /*Holds reference alingment*/
     qryAlnST,                /*Holds query alingment*/
     gapRowS,                  /*For threadsafe scoring*/
     settings                 /*Settings for alignment*/
   );

   return;
} /*HirschbergFun*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The indel column score
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      Needleman Wunsch / Smith Waterman alignment
|    o gapRowSSt to hold the last row of directions in a
|      Needleman Wunsch / Smith Waterman alignment
\--------------------------------------------------------*/
long scoreForwardHirsch(
  char *refSeqStr,        /*Reference sequence*/
  int16_t refStartS,     /*index 0 starting ref base*/
  int16_t refLenS,       /*index 1 Length of target*/

  char *qrySeqStr,        /*Query sequence*/
  int16_t qryStartS,     /*Index 0 Starting query base*/
  int16_t qryLenS,       /*index 1 length of target*/

  int16_t *scoreRowPtrS,  /*Array of scores to fill*/
  int16_t *gapRowAryS,    /*direction row, for gap extend*/
  struct alnSet *settings /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-03 TOC: scoreForwardHirsch
   '  - Does a single round of scoring for a hirschberg
   '    alignment (forward direction)
   '  o fun-03 sec-01:
   '    - Variable declerations
   '  o fun-03 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-03 sec-03:
   '    - Score till on the last row
   '  o fun-03 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *refStr = 0;
   char *refEndStr =
      refSeqStr + refStartS + refLenS - defNum16BitElms +1;

   char *qryStr = 0;
   char *qryEndStr = qrySeqStr + qryStartS + qryLenS;

   /*For sanity*/
   char bestDirC = settings->bestDirC;

   /*Loop counter variables*/
   int16_t sElm = 0;
   int16_t sStagger = 0;

   int16_t *scoreStartPtrS = scoreRowPtrS + refStartS;
   int16_t *scoreOnPtrS = scoreStartPtrS;
   int16_t insScoreS = 0;
   int16_t delScoreS = 0;
   int16_t matchScoreS = 0;

   int16_t snpBuffS[defNum16BitElms << 1]
   int16_t *snpAryS = (int16_t *) alnPointer(snpBuffS);
      
   int16_t delBuffS[defNum16BitElms << 1]
   int16_t *delAryS = (int16_t *) alnPointer(delAryS);

   int16_t nextMatchScoreL = 0;

   /*These are temporary arrays for when sequence length
   ` != vector lengths (jagged ends)*/
   int16_t tmpScoreAryS[defNum16BitElms];
   int16_t gapExtendS = gapExtendS;

   #ifndef NOGAPOPEN
      int16_t tmpGapAryS[defNum16BitElms];
      int16_t *gapStartPtrS = gapRowAryS + refStartS;
      int16_t *gapOnPtrS = gapStartPtrS;

      int16_t gapOpenS = settings->gapOpenS;
      int16_t gapDiffS = gapExtendS - gapOpenS;

      int16_t delGapBuffS[defNum16BitElms << 1]
      int16_t *delGapAryS =(int16_t *) alnPointer(delAryS);

      vectI16 gapOpenVS = set1_I16_retVectI1(gapOpenS);
      vectI16 gapDiffVS;
      vectI16 gapVS;
   #endif

   vectI16 insVS;
   vectI16 delVS;
   vectI16 snpVS;
   vectI16 gapExtendVS = set1_I16_retVectI16(gapExtendS);
   vectI16 maxVS;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-02:
   ^  - Set up the first row (indel row) of scores
   ^  o fun-03 sec-02 sub-01:
   ^    - Set up the indel column and the vector
   ^  o fun-03 sec-02 sub-02:
   ^    - Set up remaing elements (indel) (uses gap extend)
   ^  o fun-03 sec-02 sub-03:
   ^    - Fill in the last part of the indel row manually
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-01:
   *  - Set up the indel column and the vector
   \******************************************************/

   #ifndef NOGAPOPEN
      gapOnPtrS = gapRowAryS + refStartS;
      insVS = set1_I16_retVectI16(gapOpenS);
      gapVS = set1_I16_retVectI16(-1);
         /*for marking that I had gaps*/

      /*Get the difference to use in finding gap penalty*/
      gapDiffVS =
         sub_vectI16_retVectI16(gapExtendVS, gapOpenVS)
   #else
      insVS = set1_I16_retVectI16(gapExtendS);
   #endif

   indelColS = 0;

   for(
      int16_t elmInVect = 1;
      elmInVect < defNum16BitElms;
      ++elmInVect
   ){ /*Loop: build up the first part of the indel row*/
      insVS =
        add_vectI16_retVectI16(
          insVS,
          srvect_vectI16_retVectI16(gapExtendVS, elmInVect)
             /*Remove values that have been adjusted*/
      ); /*Set up the vector for filling the indel row*/
   } /*Loop: build up the first part of the indel row*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-02:
   *  - Set up indel row with vectors
   \******************************************************/

   while(
        refStr + defNum16BitElms
      < refSeqStr + refStartS + refLenS
   ){ /*Loop:Set the indel row scores (all gaps)*/
     /*Record the score*/
     storeu_vectI16_retAryI16(scoreOnPtrS, insSVS);

     #ifndef NOGAPOPEN
        /*Store that this was a gap move*/
        storeu_vectI16_retAryI16(gapOnPtrS, gapVS);
        gapOnPtrS += defNum16BitElms;
     #endif

     /*Get the next set of scores*/
     insVS = add_vectI16_retVectI16(insVS, gapExtendVS);

     scoreOnPtrS += defNum16BitElms;
     refStr += defNum16BitElms;
   } /*Loop:Set the indel row scores (all gaps)*/

   /******************************************************\
   * Fun-03 Sec-02 Sub-03:
   *  - Fill in the last part of the indel row manually
   \******************************************************/

   while(refStr < refSeqStr + refStartS + refLenS)
   { /*Loop: Fill in the end gaps in the indel row*/
      *scoreOnPtrS = *(scoreOnPtrS - 1) + gapExtendS;
      ++scoreOnPtrS;

      #ifndef NOGAPOPEN
         *gapOnPtrS = -1; /*Flag that this was a gap*/
         ++gapOnPtrS;
      #endif

      ++refStr;
   } /*Loop: Fill in the end gaps indel row*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-03:
   ^  - Set up the stagger for the vectors
   ^  o fun-04 sec-03 sub-01:
   ^    - Set starting point and run loop to find the
   ^      score for the first ref base
   ^  o fun-04 sec-03 sub-02:
   ^    - Find the scores for the first base in reference
   ^  o fun-04 sec-03 sub-03:
   ^    - Save the scores (del/snp) for the next base
   ^  o fun-04 sec-03 sub-04:
   ^    - Find the best score for the first base
   ^  o fun-04 sec-03 sub-05:
   ^    - Find the snp scores for the next base pairs
   ^  o fun-04 sec-03 sub-06:
   ^    - Find the ins scores for the next base pairs
   ^  o fun-04 sec-03 sub-07:
   ^    - Find the del scores for the next base pairs
   ^  o fun-04 sec-03 sub-08:
   ^    - Save snp/del scores needed for next base pair
   ^  o fun-04 sec-03 sub-09:
   ^    - Find the best score
   ^  o fun-04 sec-03 sub-10:
   ^    - Save snp and deletions scores for first column
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-03 Sec-03 Sub-01:
   *  - Set starting point and run loop to find the
   *    score for the first refernce base
   \******************************************************/

   scoreOnPtrS = scoreRowPtrS + refStartS;

   #ifndef NOGAPOPEN
      gapOnPtrS = gapRowAryS + refStartS;
   #endif

   /*I need to stagger the scores so that each vector
   ` is caluclating a different row
   ` For AVX512 this taks (32 is number of shorts)
   `   1 + 2 + 3 + 4 + ... + 32 = (32 * 16) rounds
   `   (16 is 32 / 2).
   */
   for(sElm = 0; sElm < defNum16BitElms; ++sElm)
   { /*Loop: stagger scores for the anti-diagnol*/
      /*First base (indel column no included)
      `  indel  1st 
      `   col   base
      ` | +0  | -10 | <--- indel row (sElm = 0) or last row
      `     \____:    Need last indel col & next indel col
      ` | -10 |<-*  |
      `               * is the cell I am finding score for
      `               + is what indelCol is set to
      ` \___ (snp), : (ins), and <- (del) show the scores
      ` needed fro each direction.
      */

      /***************************************************\
      * Fun-03 Sec-03 Sub-02:
      *  - Find the scores for the first base in reference
      \***************************************************/

      /*Find the score for an snp*/
      snpScoreS = 
          indelColS
        + getBaseScore(qrySeqStr +sElm,refSeqStr,settings);

      /*Find the next indel score*/
      #ifndef NOGAPOPEN
         /*This is only else one time, so the cpu will
         ` figure it out
         */
         if(sElm) indelColS += gapExtendS;
         else indelColS += gapOpenS;
      #else
         indelColS += gapExtendS
      #endif
     
      /*Find the deletion score*/
      delScoreS = indelColS + gapExtendS;

      /*Find the insertion score*/
      insScoreS =
          scoreOnPtrS[0]
        + (gapOnPtrS[0] & gapDiffS)
        + gapOpenS;

      /***************************************************\
      * Fun-03 Sec-03 Sub-03:
      *  - Save the scores (del/snp) for the next base
      \***************************************************/

      snpAryS[1] = delAryS[1];    /*Snp for second base*/
      delAryS[1] = scoreOnPtrS[0];/*Deletion for 2nd base*/
      delGapAryS[1] = gapOnPtrS[0];
        /*
        ` | lost|  delAryS[1]    | scoreOnPtrS[1] | *3 |
        ` | snp | scoreOnPtrS[0] |      *2        |
        ` | del |      *1        |    
        ` *1 is the cell I am filling (what scoreOnPtrS
        `   will become here (sub-02))
        ` *2 is the next base pair to score in sub-03
        ` *3 is the next next base pair to score in sub-03
        ` The idea is delAryS[1] is one base behind
        `   scoreOnPtrS. And snpAryS is two bases behind
        `   scoreOnPtrS.
        */

      /***************************************************\
      * Fun-03 Sec-03 Sub-04:
      *  - Find the best score for the first base
      \***************************************************/

      #ifndef NOGAPOPEN
         maxGapScore(
            scoreOnPtrS[0],      /*Holds best score*/
            gapOnPtrS[0],       /*Holds if was a gap*/
            snpScoreS,           /*Score for an snp*/
            insScoreS,           /*Score for an insertion*/
            insScoreS,           /*Score for an deletion*/
            settings->bestDirC   /*Prefered direction*/
         );
      #else
          alnMaxScore(
            scoreOnPtrS[0],      /*Holds best score*/
            snpScoreS,           /*Score for an snp*/
            insScoreS,           /*Score for an insertion*/
            insScoreS,           /*Score for an deletion*/
            settings->bestDirC   /*Prefered direction*/
         );
      #endif

      /***************************************************\
      * Fun-03 Sec-03 Sub-05:
      *  - Find the snp scores for the next base pairs
      \***************************************************/

      for(sStagger=1; sStagger <= sElm; ++sStagger)
      { /*Loop; Stagger each row*/
         snpScoreS += 
            getBaseScore(
               qryStr + sStagger,
               refStr + sStagger,
               settings
         ); /*Get score for the base pair*/

         snpScoreS += snpAryS[sStagger];

         /************************************************\
         * Fun-03 Sec-03 Sub-06:
         *  - Find the ins scores for the next base pairs
         \************************************************/

         #ifndef NOGAPOPEN
            insScoreS =
            ( /*find insertion score*/
                 scoreOnPtrS[sStagger]
               + (gapOnPtrS[sStagger] & gapDiffS)
               + gapOpenS
            ); /*find insertion score*/
         #else
            insScoreS = scoreOnPtrS[sStagger] + gapExtendS;
         #endif

         /************************************************\
         * Fun-03 Sec-03 Sub-07:
         *  - Find the del scores for the next base pairs
         \************************************************/

         #ifndef NOGAPOPEN
            delScoreS =
            ( /*Find the deletion score*/
                 scoreOnPtrS[sStagger - 1]
               + (gapOnPtrS[sStagger - 1] & gapDiffS)
               + gapOpenS
            ); /*Find the deletion score*/
         #else
            delScoreS =
               scoreOnPtrS[sStagger - 1] + gapExtendS;
         #endif /*Check how finding deletions*/

         /************************************************\
         * Fun-03 Sec-03 Sub-08:
         *  - Save snp/del scores needed for next base pair
         \************************************************/

         /*Score I need for the next snp*/
         snpAryS[sStagger + 1] = delAryS[sStagger + 1];
         delAryS[sStagger + 1] = scoreOnPtrS[sStagger];
         delGapAryS[sStagger + 1] = gapOnPtrS[sStagger];
         /*See sec-03 sub-03 for a diagram*/

         /************************************************\
         * Fun-03 Sec-03 Sub-09:
         *  - Find the best score
         \************************************************/

         #ifndef NOGAPOPEN
            maxGapScore(
               scoreOnPtrS[sStagger],/*Holds best score*/
               gapOnPtrS[sStagger],/*Holds if was a gap*/
               snpScoreS,        /*Score for an snp*/
               insScoreS,        /*Score for an insertion*/
               delScoreS,        /*Score for an deletion*/
               settings->bestDirC/*Prefered direction*/
            );
         #else
             alnMaxScore(
               scoreOnPtrS[sStagger],/*Holds best score*/
               snpScoreS,        /*Score for an snp*/
               insScoreS,        /*Score for an insertion*/
               insScoreS,        /*Score for an deletion*/
               settings->bestDirC/*Prefered direction*/
            );
         #endif
      } /*Loop; Stagger each row*/
   } /*Loop: stagger scores for the anti-diagnol*/

   /*******(**********************************************\
   * Fun-03 Sec-03 Sub-10:
   *  - Save snp and deletions scores for first column
   \******************************************************/

   snpAryS[0] = indelColS - gapExtendS;
   delAryS[0] = indelColS;
   delGapAryS[0] = -1; /*This is always a gap (indel col)*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-04:
   ^  - Find the scores for each row until at the end
   ^  o fun-03 sec-04 diag-01:
   ^    - diagram of staggard matrix before using vectors
   ^  o fun-03 sec-04 sub-01:
   ^    - Load up intitial scores into vectors &start loops
   ^  o fun-03 sec-04 sub-02:
   ^    - Find the snp/match scores for this round
   ^  o fun-03 sec-04 sub-03:
   ^    - Find the deletion scores for this round
   ^  o fun-03 sec-04 sub-04:
   ^    - Add in new bases insertion score for this round
   ^  o fun-03 sec-04 sub-05:
   ^    - Add in if the new base is an gap extend or open
   ^  o fun-03 sec-04 sub-06:
   ^    - Select best score for each base pair
   ^  o fun-03 sec-04 diag-02
   ^    - diagram of matrix when starting a new row
   ^  o fun-03 sec-04 sub-07:
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /* Fun-03 Sec-04 Diag-01:
   `  What the matrix looks like * = found, o = missing,
   `  g = gap column/row, ... or (.) = repeat (Only first
   `  conor filled)
   ` |g|g|g|g|g|g|g|g|g|g|g|g|...|g|
   ` |g|*|*|*|*|*|*|*|*|o|o|o|...|o|
   ` |g|*|*|*|*|*|*|*|o|o|o|o|...|o|
   ` |g|*|*|*|*|*|*|o|o|o|o|o|...|o|
   ` |g|*|*|*|*|*|o|o|o|o|o|o|...|o|
   ` |g|*|*|*|*|o|o|o|o|o|o|o|...|o|
   ` |g|*|*|*|o|o|o|o|o|o|o|o|...|o|
   ` |g|*|*|o|o|o|o|o|o|o|o|o|...|o|
   ` |g|*|o|o|o|o|o|o|o|o|o|o|...|o|
   ` |g|o|o|o|o|o|o|o|o|o|o|o|...|o|
   `  . . . . . . . . . . . . ... .
   `  . . . . . . . . . . . . ... .
   `  . . . . . . . . . . . . ... .
   ` |g|o|o|o|o|o|o|o|o|o|o|o|...|o|
   */   

   /******************************************************\
   * Fun-04 Sec-04 Sub-01:
   *  - Load up intitial scores into vectors & start loops
   \******************************************************/

   /*This the setup the loop will use*/
   maxVS = loadu_vectI16(scoreOnPtrS); /*Has found scores*/
   insVS = load_vectI16_retVectI16(snpAryS);

   #ifndef NOGAPOPEN
      gapVS = loadu_I16Ary_retVectI16(gapOnPtrS);
   #endif

  /*The qryEndStr and refEndStr already have
   `  defNum16BitElms accounted into them
   */

   while(qryStr < qryEndStr)
   { /*Loop: score all query bases (rows)*/
     refStr = refSeqStr + refStartS;

     while(refStr < refEndStr)
     { /*Loop:score all query bases (columns)*/

       /**************************************************\
       * Fun-04 Sec-04 Sub-02:
       *  - Find the snp/match scores for this round
       \**************************************************/

       /*I need to find the match/mismatch scores*/
       for(sElm = 0; sElm < defNum16Elms; ++sElm)
       { /*Loop: Get scores for if matches or not*/
         snpAryS[sElm] =
          getBaseScore(qryStr[sElm],refStr[sElm],settings);
       } /*Loop: Get scores for if matches or not*/

       ++refStr;

       /*Using delVS as a temporary vector*/
       delVS = load_aryI16_retVectI16(snpAryS);
       snpScoreVS = add_vectI16_retVectI16(delVS, insVS);

       /**************************************************\
       * Fun-04 Sec-04 Sub-03:
       *  - Find the deletion scores for this round
       \**************************************************/

        /*last rounds max scores = this rounds deletions*/

        #ifndef NOGAPOPEN
          delVS = and_vectI16_retVectI16(gapVS, gapDiffVS);
          delVS = add_vectI16_retVectI16(delVS, gapOpenVS);
          delVS = add_vectI16_retVectI16(maxVS, delVS);
        #else
          delVS=add_vectI16_retVectI16(maxVS,gapExtendVS);
        #endif

       /**************************************************\
       * Fun-04 Sec-04 Sub-04:
       *  - Add in new bases insertion score for this round
       \**************************************************/

        /*Insertion*/
        *scoreOnPtr =
           extract_vectI16_I16(maxVS, defNum16BitElms - 1);
           /*Vectors load backwards, so the element that I
           ` discarding is at the end of the vector, not
           ` the start
           */

        maxVS = slvect_vectI16_retVectI16(maxVS, 1);
           /*Vectors load backwards, so this should be
           ` shifting the frist score out
           */
        maxVS =
           insert_I16_retVectI16(
              maxVS,
              *(scoreOnPtr + defNum16BitElms - 1),
              0 /*Vectors load backwards,so last is first*/
        ); /*Insert the new score*/

        ++scoreOnPtr; /*Next score to store next round*/

       /**************************************************\
       * Fun-04 Sec-04 Sub-05:
       *  - Add in if the new base is an gap extend or open
       \**************************************************/

        #ifndef NOGAPOPEN
           *gapOnPtrS =
              extract_vectI16_I16(gapVS,defNum16BitElms-1);

           gapVS = slvect_vectI16_retVectI16(gapVS, 1);

           gapVS =
              insert_I16_retVectI16(
                 gapVS,
                 *(gapOnPtrS + defNum16BitElms - 1),
                 0 /*Vectors load backwards,so last is first*/
           ); /*Insert the new gap*/

           ++gapOnPtrS; /*Next gap to store next round*/
        #endif

       /**************************************************\
       * Fun-04 Sec-04 Sub-05:
       *  - Find the insertion score for the new base
       \**************************************************/

        #ifndef NOGAPOPEN
          insVS = and_vectI16_retVectI16(gapVS, gapDiffVS);
          insVS = add_vectI16_retVectI16(insVS, gapOpenVS);
          insVS = add_vectI16_retVectI16(maxVS, insVS);
        #else
          insVS=add_vectI16_retVectI16(insVS,gapExtendVS);
        #endif

       /**************************************************\
       * Fun-03 Sec-04 Sub-06:
       *  - Select best score for each base pair
       \**************************************************/

       #if !defined NOGAPOPEN
        vI16MaxGap(maxVS,gapVS,snpVS,insVS,delVS,bestDirC);
       #else
        vI16Max(maxVS, snpVS, insVS, delVS, bestDirC);
       #endif
     } /*Loop:score all query bases (columns)*/

     /* Fun-03 Sec-04 Diag-02:
     ` At this point I am at the end of the row and need
     ` to handle starting a new row each round. The first
     ` matrix shows the state I am at. The second matrix
     ` showws the state after I find the next maximum.
     ` * = max score is found, o = need fo find max score.
     `  |*|...|*|*|*|*|*|*|*|*|    |*|...|*|*|*|*|*|*|*|*|
     `  |*|...|*|*|*|*|*|*|*|o|    |*|...|*|*|*|*|*|*|*|*|
     `  |*|...|*|*|*|*|*|*|o|o|    |*|...|*|*|*|*|*|*|*|o|
     `  |*|...|*|*|*|*|*|o|o|o|    |*|...|*|*|*|*|*|*|o|o|
     `  |*|...|*|*|*|*|o|o|o|o|--->|*|...|*|*|*|*|*|o|o|o|
     `  |*|...|*|*|*|o|o|o|o|o|    |*|...|*|*|*|*|o|o|o|o|
     `  |*|...|*|*|o|o|o|o|o|o|    |*|...|*|*|*|o|o|o|o|o|
     `  |*|...|*|o|o|o|o|o|o|o|    |*|...|*|*|o|o|o|o|o|o|
     `  |*|...|o|o|o|o|o|o|o|o|    |*|...|*|o|o|o|o|o|o|o|
     `  |o|...|o|o|o|o|o|o|o|o|    |*|...|o|o|o|o|o|o|o|o|
     `  |o|...|o|o|o|o|o|o|o|o|    |o|...|o|o|o|o|o|o|o|o|
     ` I need to account for startring a new row for each
     ` element in the vector (8 for 128 bit vector).
     `
     ` The first score has snp with the previous indel
     ` column. The deletion is the for the new score is
     ` new indel column (not shown in diagram)
     */

     sStagger = 0;

     while(refStr < refSeqStr + refStartS + refLenS)
     { /*Loop: Find the last scores in the row*/

       /**************************************************\
       * Fun-04 Sec-04 Sub-07:
       *  - Find the snp/match scores for the rows end
       \**************************************************/

       /*I need to find the match/mismatch scores*/
       for(sElm = 0; sElm < sStagger; ++sElm)
       { /*Loop: Get scores for if matches or not*/
         snpAryS[sElm] =
          getBaseScore(qryStr[sElm],refStr[sElm],settings);
       } /*Loop: Get scores for if matches or not*/

       snpAryS[sElm] =
          getBaseScore(
             qryStr[defNum16BitElms],
             refStr[sElm],
             settings
       );

       for(sElm = sStagger+1; sElm < defNum16Elms; ++sElm)
       { /*Loop: Get scores for if matches or not*/
         snpAryS[sElm] =
          getBaseScore(qryStr[sElm],refStr[sElm],settings);
       } /*Loop: Get scores for if matches or not*/

       /*The first ref base always has indelCol as snp*/
       snpAryS[sStagger] += indelColS;
       indelColS += gapExtendS; /*Score for next row*/

       ++refStr;

       /*Make sure the new snp score is new. sStagger
       ` is used without correction because vectors are
       ` loaded backwards
       */
       insVS = insert_I16_retVectI16(insVS, 0, sStagger);

       /*Using delVS as a temporary vector*/
       delVS = load_aryI16_retVectI16(snpAryS);
       snpScoreVS = add_vectI16_retVectI16(delVS, insVS);

       /**************************************************\
       * Fun-04 Sec-04 Sub-08:
       *  - Find the deletion scores for the rows end
       \**************************************************/

       /*I need to remove the row end score/gap*/
       *(scoreOnPtrS + defNum16BitElms - 1) =
          extract_vectI16_retI16(maxVS,sStagger);

       maxVS =
          insert_I16_retVectI16(maxVS,indelColS,sStagger);

       #ifndef NOGAPOPEN
         *(gapOnPtrS - defNum16BitElms) =
            extract_vectI16_retI16(gapVS, sStagger);

         /*This is always a gap extension penalty*/
         gapVS = insert_I16_retVectI16(gapVS,-1,sStagger);

         delVS = and_vectI16_retVectI16(gapVS, gapDiffVS);
         delVS = add_vectI16_retVectI16(delVS, gapOpenVS);
         delVS = add_vectI16_retVectI16(maxVS, delVS);
       #else
         delVS=add_vectI16_retVectI16(maxVS,gapExtendVS);
       #endif

       /**************************************************\
       * Fun-04 Sec-04 Sub-09:
       *  - Find the insertion scores for the rows end
       \**************************************************/

       maxVS =
          insert_I16_retVectI16(
             maxVS,
             *(scoreStartPtrS + sStagger),
             sStagger
       ); /*Insert the new rows insertion score*/

       #ifndef NOGAPOPEN
          gapVS =
             insert_I16_retVectI16(
                gapVS,
                *(gapStartPtrS + sStagger),
                sStagger
          ); /*Insert if the insertion is a gap extend*/
       #endif

       /**************************************************\
       * Fun-04 Sec-04 Sub-10:
       *  - Add in new bases insertion score for row end
       \**************************************************/

        *scoreOnPtr =
           extract_vectI16_I16(maxVS, defNum16BitElms - 1);
           /*Vectors load backwards, so the element that I
           ` discarding is at the end of the vector, not
           ` the start
           */

        maxVS = slvect_vectI16_retVectI16(maxVS, 1);
           /*Vectors load backwards, so this should be
           ` shifting the frist score out
           */
        maxVS =
           insert_I16_retVectI16(
              maxVS,
              *(scoreOnPtr + defNum16BitElms - 1),
              0 /*Vectors load backwards,so last is first*/
        ); /*Insert the new score*/

        ++scoreOnPtr; /*Next score to store next round*/

       /**************************************************\
       * Fun-04 Sec-04 Sub-11:
       *  - Add in if base is gap extend or open row end
       \**************************************************/

        #ifndef NOGAPOPEN
           *gapOnPtrS =
              extract_vectI16_I16(gapVS,defNum16BitElms-1);

           gapVS = slvect_vectI16_retVectI16(gapVS, 1);

           gapVS =
              insert_I16_retVectI16(
                 gapVS,
                 *(gapOnPtrS + defNum16BitElms - 1),
                 0/*Vectors load backwards,so last is 1st*/
           ); /*Insert the new gap*/

           ++gapOnPtrS; /*Next gap to store next round*/
        #endif

       /**************************************************\
       * Fun-04 Sec-04 Sub-12:
       *  - Find insertion score for new base near row end
       \**************************************************/

        #ifndef NOGAPOPEN
          insVS = and_vectI16_retVectI16(gapVS, gapDiffVS);
          insVS = add_vectI16_retVectI16(insVS, gapOpenVS);
          insVS = add_vectI16_retVectI16(maxVS, insVS);
        #else
          insVS=add_vectI16_retVectI16(insVS,gapExtendVS);
        #endif

       /**************************************************\
       * Fun-03 Sec-04 Sub-13:
       *  - Select best score for each base pair row end
       \**************************************************/

       #if !defined NOGAPOPEN
        vI16MaxGap(maxVS,gapVS,snpVS,insVS,delVS,bestDirC);
       #else
        vI16Max(maxVS, snpVS, insVS, delVS, bestDirC);
       #endif

       ++sStagger;
     } /*Loop: Find the last scores in the row*/

     /****************************************************\
     * Fun-03 Sec-04 Sub-14:
     *  - Prepare for next loop
     \****************************************************/

     scoreOnPtrL = scoreStartPtrS;

     #ifndef NOGAPOPEN
        gapOnPtr = gapStartPtrS;
     #endif

     qryStr += defNum16BitElms;
   } /*Loop: score all query bases (rows)*/
  
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-03 Sec-05:
   ^  - Unstagger (find) the last scores in the conor
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
 
   /* What the matrix looks like * = found, o = missing,
   `  ... or (.) = repeat (Only last conor not filled)
   ` |*|*|*|...|*|*|*|*|*|*|*|*|*|
   ` |*|*|*|...|*|*|*|*|*|*|*|*|*|
   ` |*|*|*|...|*|*|*|*|*|*|*|*|*|
   `  . . . ... . . . . . . . . .
   `  . . . ... . . . . . . . . .
   `  . . . ... . . . . . . . . .
   ` |*|*|*|...|*|*|*|*|*|*|*|*|o|
   ` |*|*|*|...|*|*|*|*|*|*|*|o|o|
   ` |*|*|*|...|*|*|*|*|*|*|o|o|o|
   ` |*|*|*|...|*|*|*|*|o|o|o|o|o|
   ` |*|*|*|...|*|*|*|o|o|o|o|o|o|
   ` |*|*|*|...|*|*|o|o|o|o|o|o|o|
   ` |*|*|*|...|*|o|o|o|o|o|o|o|o|
   */   


  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-03 Sec-04:
  ^  - Clean up
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Correct for being on the next row*/
   indelColL -= settings->gapExtendS;
   return indelColL;
} /*scoreForwardHirsch*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o The indel column score
|  - Modifies:
|    o scoreRowPtrL to hold the last row of scores in a
|      backwards Needleman Wunsch /Smith Waterman alignment
|    o gapRowSST to hold the last row of directions in a
|      backwards Needleman Wunsch /Smith Waterman alignment
\--------------------------------------------------------*/
long scoreReverseHirsch(
  char *refSeqStr,          /*Reference sequence*/
  unsigned long refStartUL,  /*index 0 starting ref base*/
  unsigned long refLenUL,    /*index 1 Length of target*/

  char *qrySeqStr,          /*Query sequence*/
  unsigned long qryStartUL, /*Index 0 Starting query base*/
  unsigned long qryLenUL,    /*index 1 length of target*/

  long *scoreRowPtrL,        /*Array of scores to fill*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *gapRowSST,/*direction row*/
  #else
     char *gapRowSST,      /*direction row, for gap extend*/
  #endif
  struct alnSet *settings    /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-04 TOC: scoreReverseHirsch
   '  - Does a single round of scoring for a hirschberg
   '    alignment (reverse direction)
   '  o fun-04 sec-01:
   '    - Variable declerations
   '  o fun-04 sec-02:
   '    - Set up the first row (indel row) of scores
   '  o fun-04 sec-03:
   '    - Score till on the last row
   '  o fun-04 sec-04:
   '    - Clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   char *qryStr = 0;
   char *refStr = 0;

   long *scoreOnPtrL = 0;
     /*refLenUL is index 1, so  need -1 to get index 0*/

   long insScoreL = 0;
   long delScoreL = 0;
   long matchScoreL = 0;
   long nextMatchScoreL = 0; /*Score for the next base*/
   long indelColL = 0;       /*Holds indel column values*/

   #if !defined HIRSCHTWOBIT && !defined NOGAPOPEN
      char *endDirRowStr = 0;
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-02:
   ^  - Set up the first row (indel row) of scores
   ^  o fun-04 sec-02 sub-01:
   ^    - Set up the first two elements (no gat extend)
   ^  o fun-04 sec-02 sub-02:
   ^    - Set up remaing elements (indel) (uses gap extend)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-02 Sub-01:
   *  - Set up the first two elements (no gat extend)
   \******************************************************/

   #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
     twoBitMvXElmFromStart(gapRowSST,refStartUL+refLenUL-1);
     changeTwoBitElm(gapRowSST, defMvDel);
     twoBitMvBackOneElm(gapRowSST);
   #elif !defined NOGAPOPEN
     gapRowSST += refStartUL + refLenUL - 1;
     endDirRowStr = gapRowSST;
     *gapRowSST = defMvDel;
     --gapRowSST;
   #endif

   scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL - 1;
     /*- 1 to account for refLenUL being index 1*/
   indelColL = 0;

   #if !defined NOGAPOPEN
      *scoreOnPtrL = settings->gapOpenS;
   #else
      *scoreOnPtrL = settings->gapExtendS;
   #endif

   --scoreOnPtrL;

   /******************************************************\
   * Fun-04 Sec-02 Sub-02:
   *  - Set up remaing elements (indels) (uses gap extend)
   \******************************************************/

   /* Loop from the second to last base till the start,
   `  Starting at second to last, because already did the
   `  first base
   */
   refStr = refSeqStr + refStartUL + refLenUL - 2;

   while(refStr > refSeqStr + refStartUL - 1)
   { /*Loop:Set the initial blank scores*/
     *scoreOnPtrL = *(scoreOnPtrL+1) +settings->gapExtendS;

     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        changeTwoBitElm(gapRowSST, defMvDel);
        twoBitMvBackOneElm(gapRowSST);
     #elif !defined NOGAPOPEN
        *gapRowSST = defMvDel;
        --gapRowSST;
     #endif

     --scoreOnPtrL;
     --refStr;
   } /*Loop:Set the initial blank scores*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-04 Sec-03:
   ^  - Score till on the last row
   ^  o fun-04 sec-03 sub-01:
   ^    - Set up the first element (indel, has gap open)
   ^  o fun-04 sec-03 sub-02:
   ^    - Find the scores for the next row
   ^  o fun-04 sec-03 sub-03:
   ^    - Select best score for direction
   ^  o fun-04 sec-03 sub-04:
   ^    - Find the scores for the next indels
   ^  o fun-04 sec-03 sub-05:
   ^    - Find the best score for last base pair in row
   ^  o fun-04 sec-03 sub-06:
   ^    - Move to start of row (direction matrix)
   ^  o fun-04 sec-03 sub-07:
   ^    - Find the scores for the 1st baise pair in the row
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-04 Sec-03 Sub-01:
   *  - Set up the first element (indel, has gap open)
   \******************************************************/

   /*Move back to the start of row (refLenUL is index 1)*/
   #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
     twoBitMvXElmFromStart(gapRowSST,refStartUL+refLenUL-1);
   #elif !defined NOGAPOPEN
     gapRowSST = endDirRowStr;
   #endif

   scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL - 1;

   /*Find the first match/snp score (first ref base)*/
   nextMatchScoreL =
     getBaseScore(
       qrySeqStr + qryStartUL + qryLenUL - 1,
       refSeqStr + refStartUL + refLenUL - 1,
       settings
   ); /*Get the score for the next base*/

   nextMatchScoreL += indelColL;
     /* This is here so I can overwrite the array with the
     `  new scores
     */

   /*Find the first insertion and deletion scores*/

   #if !defined NOGAPOPEN
      indelColL += settings->gapOpenS;
   #else
      indelColL += settings->gapExtendS;
   #endif

   insScoreL = indelColL + settings->gapExtendS;
   delScoreL = indelColL + settings->gapExtendS;

   /******************************************************\
   * Fun-04 Sec-03 Sub-02:
   ^  - Find the scores for the next row
   \******************************************************/

   refStr = refSeqStr + refStartUL + refLenUL - 2;
     /*Start one base back to account for the
     ` pre-calcualted scores for the first base
     */
   qryStr = qrySeqStr + qryStartUL + qryLenUL - 1;

   while(qryStr > qrySeqStr + qryStartUL - 1)
   { /*Loop: score all query bases (rows)*/

     while(refStr > refSeqStr + refStartUL - 1)
     { /*Loop:score all query bases (columns)*/

       /* Get the score for the next match. This allows me
       `  to overwite the last diagnol and thus, use only
       `  one row of scoring
       */

       matchScoreL = nextMatchScoreL;

       nextMatchScoreL =
         getBaseScore(
           qryStr,     /*Current query base*/
           refStr,     /*next reference base*/
           settings    /*Has score matrix*/
       ); /*Get the score for the enxt base*/
          /*At worst case refStr will be '\0'*/

       nextMatchScoreL += *scoreOnPtrL;

       /**************************************************\
       * Fun-04 Sec-03 Sub-03:
       *  - Select best score for direction
       \**************************************************/

       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          twoBitMaxScore(
            gapRowSST,
            settings,
            &insScoreL,
            &matchScoreL,
            &delScoreL,
            scoreOnPtrL
          ); /*Update the score and direction*/
       #elif !defined NOGAPOPEN
          charMaxScore(
            gapRowSST,        /*Direction matrix*/
            settings,        /*has direction preference*/
            &insScoreL,      /*Score for an insertion*/
            &matchScoreL,    /*Score for an deletion*/
            &delScoreL,      /*Score for an match/snp*/
            scoreOnPtrL      /*Score position to update*/
          ); /*Update the score and direction*/
       #else
          alnMaxScore(
             settings,
             &insScoreL,
             &matchScoreL,
             &delScoreL,
             scoreOnPtrL
          );
       #endif

       /**************************************************\
       * Fun-04 Sec-03 Sub-04:
       *  - Find the scores for the next indels
       \**************************************************/

       /* The deletion scores are based on the found base,
       `  So I can find the next score before moving
       `  Get the deletion score
       */
       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          indelScore(
             delScoreL,
             getTwoBitElm(gapRowSST),
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/

          twoBitMvBackOneElm(gapRowSST);
       #elif !defined NOGAPOPEN
         indelScore(
             delScoreL,
             *gapRowSST, /*For no gap open, is a dummy*/
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/

          --gapRowSST;
       #else
         delScoreL = *scoreOnPtrL + settings->gapExtendS;
      #endif

       --scoreOnPtrL;

       /* Finding indel scores at end, so that I can keep
       `  the indel column in a separate variable
       `  Get the insertion score
       */
       #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
          indelScore(
             insScoreL,
             getTwoBitElm(gapRowSST),
             *scoreOnPtrL,
             settings
          ); /*Macro from generalAlnFun.h*/
       #elif !defined NOGAPOPEN
          indelScore(
              insScoreL,
              *gapRowSST, /*For no gap open, is a dummy*/
              *scoreOnPtrL,
              settings
           ); /*Macro from generalAlnFun.h*/
       #else
          insScoreL = *scoreOnPtrL + settings->gapExtendS;
       #endif

       --refStr;
     } /*Loop:score all query bases (columns)*/

     /****************************************************\
     * Fun-04 Sec-03 Sub-05:
     *  - Find the best score for last base pair in row
     \****************************************************/

     /*Update the final score in the row*/
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        twoBitMaxScore(
          gapRowSST,
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        ); /*Update the score and direction*/
     #elif !defined NOGAPOPEN
        charMaxScore(
          gapRowSST,
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        ); /*Update the score and direction*/
     #else
       alnMaxScore(
          settings,
          &insScoreL,
          &nextMatchScoreL,
          &delScoreL,
          scoreOnPtrL
        );
     #endif

     /****************************************************\
     * Fun-04 Sec-03 Sub-06:
     *  - Move to start of row (direction matrix)
     \****************************************************/

     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
       twoBitMvXElmFromStart(
          gapRowSST,
          refStartUL + refLenUL - 1
       );
     #elif !defined NOGAPOPEN
        gapRowSST = endDirRowStr;
     #endif

     scoreOnPtrL = scoreRowPtrL + refStartUL + refLenUL -1;

     /****************************************************\
     * Fun-04 Sec-03 Sub-07:
     *  - Find the scores for the 1st baise pair in the row
     \****************************************************/

     /*I need to refind the insertion and deletion scores*/
     #if defined HIRSCHTWOBIT && !defined NOGAPOPEN
        indelScore(
           insScoreL,
           getTwoBitElm(gapRowSST),
           *scoreOnPtrL,
           settings
        ); /*Macro from generalAlnFun.h*/
     #elif !defined NOGAPOPEN
        indelScore(
            insScoreL,
            *gapRowSST, /*For no gap open, is a dummy*/
            *scoreOnPtrL,
            settings
         ); /*Macro from generalAlnFun.h*/
     #else
        insScoreL = *scoreOnPtrL + settings->gapExtendS;
     #endif

     /*Reset the base on*/
     refStr = refSeqStr + refStartUL + refLenUL - 2;
     --qryStr;

     /*Find the first match/snp score (first ref base)*/
     nextMatchScoreL =
       getBaseScore(
         qryStr,             /*Next query base*/
         refStr + 1,
         settings            /*Has score matrix*/
     ); /*Get the score for the enxt base*/

     nextMatchScoreL += indelColL;
     indelColL += settings->gapExtendS;
     delScoreL = indelColL+settings->gapExtendS;
   } /*Loop: score all query bases (rows)*/
  
  /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
  ^ Fun-04 Sec-04:
  ^  - Clean up
  \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*Correct for being on the next row*/
   indelColL -= settings->gapExtendS;
   return indelColL;
} /*scoreReverseHirsch*/

/*--------------------------------------------------------\
| Output:
|  - Modifies:
|    o twoBitAlnST to hold the alignment for the single
|      base aligned to the sequence
\--------------------------------------------------------*/
void vectI16PosOneBase(
  char baseC,           /*Single base to align to a seq*/
  int16_t baseIndexS,   /*position of baseC; Index 0*/
  char *seqStr,         /*Sequence position on*/
  int16_t startOfSeqS,  /*index 0; 1st base to align*/
  int16_t lenSeqS,      /*seqStr target length (index 1)*/
  char *baseCAlnST,     /*Alignemt array for baseC*/
  char *seqAlnST,       /*Alignment array for seqStr*/
  struct alnSet *settings /*setttings to use*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-05 TOC: positionSingleRefBase
   '  - Align a single base to a sequence
   '  o fun-05 sec-01:
   '    - Variable declerations
   '  o fun-05 sec-02:
   '    - Find the reference bases position on the query
   '  o fun-05 sec-03:
   '    - Fill in insertions and reference base position
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   int16_t insScoreS = 0;    /*The score of an insertion*/
   int16_t delScoreS = 0;    /*the score of a deleton*/
   int16_t matchScoreS = 0;  /*The score of a match*/
   int16_t curScoreS = 0;    /*The score of the last row*/

   char *seqBaseStr = seqStr + startOfSeqS;
   char *endSeqStr = seqStr +startOfSeqS +lenSeqS-1;
   char *lastMatchStr = 0;
   char dirC = 0; /*Removes some #if defined statments*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-02:
   ^  - Find the reference bases position on the query
   ^  o fun-05 sec-02 sub-01:
   ^    - Find the first scores for the loop
   ^  o fun-05 sec-02 sub-02:
   ^    - Find the remaing scores
   ^  o fun-05 sec-02 sub-03:
   ^     - Figure out which of hte final scores to keep
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-05 Sec-02 Sub-01:
   *  - Find the first scores for the loop
   \******************************************************/

   baseCAlnST += baseIndexS;
   seqAlnST += startOfSeqS;
   matchScoreS = getBaseScore(seqBaseStr,&baseC,settings);

   #if !defined NOGAPOPEN
      insScoreS = settings->gapOpenS;
   #else
      insScoreS = settings->gapExtendS;
   #endif

   delScoreS = insScoreS;
   ++seqBaseStr;

   /******************************************************\
   * Fun-05 Sec-02 Sub-02:
   *  - Find the remaing scores
   \******************************************************/
 
   do { /*While I have bases to compare*/
     charMaxScore(
       &dirC,
       settings,        /*has direction preference*/
       &insScoreS,      /*Score for an insertion*/
       &matchScoreS,    /*Score for an deletion*/
       &delScoreS,      /*Score for an match/snp*/
       &curScoreS       /*Score position to update*/
     ); /*Update the score*/
     /*Directionalty determines which direction to select
     ` when one or more directions are equal (ins, snp,del)
     */

     matchScoreS =
         insScoreS
       + getBaseScore(seqBaseStr,&baseC,settings);

     insScoreS += (int16_t) settings->gapExtendS;

     switch(dirC)
     { /*Switch: Check if keeping the score*/
       case defMvSnp:
         lastMatchStr = seqBaseStr - 1;
         #if !defined NOGAPOPEN
            delScoreS =
               curScoreS + (int16_t) settings->gapOpenS;
         #else
            delScoreS =
               curScoreS + (int16_t) settings->gapExtendS;
         #endif
         break;

       case defMvDel:
       case defMvIns:
         delScoreS =
            curScoreS + (int16_t) settings->gapExtendS;
         break;
       case defMvStop: break; /*Never fires*/
     } /*Switch: Check if keeping the score*/

     ++seqBaseStr;
   } while(seqBaseStr <= endSeqStr);
   /*While I have bases to compare*/

   /******************************************************\
   * Fun-05 Sec-02 Sub-03:
   *  - Figure out which of the final scores to keep
   \******************************************************/

   charMaxScore(
     &dirC,
     settings,        /*has direction preference*/
     &insScoreS,      /*Score for an insertion*/
     &matchScoreS,    /*Score for an deletion*/
     &delScoreS,      /*Score for an match/snp*/
     &curScoreS       /*Score position to update*/
   ); /*Update the final score*/
   /*Directionalty determines which direction to select
   ` when one or more directions are equal (ins, snp,del)
   */

   switch(dirC)
   { /*Switch: Check if keeping the score*/
     case defMvSnp:
       lastMatchStr = seqBaseStr - 1;
       break;

     case defMvDel: break;
     case defMvIns: break;
     case defMvStop: break; // Never fires
   } /*Switch: Check if keeping the score*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-05 Sec-03:
   ^  - Fill in the insertions and reference base position
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   /* A series of deletions and insertions are prefered
   `  over matches and smps. In this case put the base
   `  at the first base. There is no good position
   */
   if(lastMatchStr == 0)
     lastMatchStr = seqStr + startOfSeqS;

   seqBaseStr = seqStr + startOfSeqS;
   /* No need to change query position since previouis
   `  loop only usedo one direction position
   */

   vectAddGaps(
      seqAlnST,
      seqBaseStr - lastMatchStr, /*Index 1*/
      defGapFlag
   ); /*Add in gaps before the last match*/

   seqBaseStr = lastMatchStr;
   
   /*Add in the position of the base*/
   if(checkIfBasesMatch(seqBaseStr, &baseC))
   { /*IF the bases matched*/
      *baseCAlnST = defMatchFlag;
      *seqAlnST = defMatchFlag;
   } /*IF the bases matched*/

   else
   { /* Else this was a SNP*/
      *baseCAlnST = defSnpFlag;
      *seqAlnST = defSnpFlag;
   } /*Else this was a SNPS*/

   ++seqAlnST;
   ++seqBaseStr;

   vectAddGaps(
      seqAlnST,
      1 + endSeqStr - seqBaseStr, /*+1 to make index 1*/
      defGapFlag
   ); /*Add in gaps after the last match*/

   return;
} /*positionSingleRefBase*/

/*--------------------------------------------------------\
| Output:
|  - Returns:
|    o 0: for error
|    o pointer to alnStruct with alignment
\--------------------------------------------------------*/
struct alnStruct * hirschVectToAlnST(
  struct seqStruct *refST,
   /*Has reference alignment start/end & reference length*/
  struct seqStruct *qryST,
   /*Has query alignment start/end and query length*/
  #ifdef HIRSCHTWOBIT
     struct twoBitAry *refAlignment,
       /*Two bit array with the reference alignment*/
     struct twoBitAry *qryAlignment
       /*Two bit array with the query alignment*/
  #else
     char *refAlignment, /*has reference alignment*/
     char *qryAlignment /*has query alignment*/
  #endif
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-06 TOC: twoBitAlnToAlnST
   '  - Converts a two bit array with an alignment to an
   '    alnStruct structure
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   char *refAlnStr = 0;
   char *qryAlnStr = 0;
   uint8_t bitUC = 0;

   long refIndexUL = 0;
   long qryIndexUL = 0;

   long refFirstAlnBaseL = -1;
   long refLastAlnBaseL = 0;
   long qryFirstAlnBaseL = -1;
   long qryLastAlnBaseL = 0;

   long numDelsL = 0;
   long numInssL = 0;
   long numSnpsL = 0;
   long numMatchesL = 0;

   struct alnStruct *alnST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-02:
   ^  - Allocate memory
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   alnST = malloc(sizeof(alnStruct));
   if(alnST == 0) return 0; /*Memory error*/
   initAlnST(alnST);

   refAlnStr = calloc(refST->lenSeqUL + 1, sizeof(char));

   if(refAlnStr == 0)
   { /*If had a memor error*/
      freeAlnST(alnST, 1); /*1 for freeing heap memory*/
      alnST = 0;
      return 0;
   } /*If had a memor error*/

   alnST->refAlnStr = refAlnStr;

   qryAlnStr = calloc(qryST->lenSeqUL + 1, sizeof(char));

   if(qryAlnStr == 0)
   { // IF had a memory allocation error
     freeAlnST(alnST, 1); /*1 for freeing heap memory*/
     alnST = 0;
     refAlnStr = 0;
     return 0;
   } // IF had a memory allocation error

   alnST->qryAlnStr = qryAlnStr;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-03:
   ^  - Add in the alignment
   ^  o fun-06 sec-03 sub-01:
   ^    - Add in the reference alignment
   ^  o fun-06 sec-03 sub-02:
   ^    - Add in the query alignment
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Add in the reference alignment
   \******************************************************/

   #ifdef HIRSCHTWOBIT
      twoBitMvXElmFromStart(refAlignment, 0);
   #endif

   /*the refAlignment array starts at refST->offsetUL*/
   refIndexUL = refST->offsetUL;
   refAlnStr += refIndexUL;

   #ifdef HIRSCHTWOBIT
      bitUC = getTwoBitElm(refAlignment);
   #else
      bitUC = (uint8_t) *refAlignment;
   #endif

   while(bitUC != defEndAlnFlag)
   { /*Loop: Add reference aligned bases to alnStruct*/
      switch(bitUC)
      { /*Switch: Check the error type*/
         case 0: break;

         case defSnpFlag:
         /*Case: snp*/
            ++numSnpsL;

            if(refFirstAlnBaseL < 0) /*Start of alignment*/
               refFirstAlnBaseL = refIndexUL;

            refLastAlnBaseL=refIndexUL;/*end of alignment*/
            break;
         /*Case: snp*/

         case defMatchFlag:
         /*Case: matches*/
            ++numMatchesL;

            if(refFirstAlnBaseL < 0) /*Start of alignment*/
               refFirstAlnBaseL = refIndexUL;

            refLastAlnBaseL=refIndexUL;/*end of alignment*/
            break;
         /*Case: matches*/

         case defGapFlag:
         /*Case: Deletions*/
            ++numDelsL;
            break;
         /*Case: Deletions*/
      } /*Switch: Check the error type*/

      *refAlnStr = bitUC;

      #ifdef HIRSCHTWOBIT
         twoBitMvToNextElm(refAlignment);
         bitUC = getTwoBitElm(refAlignment);
      #else
         ++refAlignment;
         bitUC = (uint8_t) *refAlignment;
      #endif

      ++refAlnStr;
      ++refIndexUL;
   } /*Loop: Add reference aligned bases to alnStruct*/

   if(refFirstAlnBaseL >= 0) /*Start of alignment*/
      alnST->refStartAlnUL = refFirstAlnBaseL;
   else alnST->refStartAlnUL = refST->lenSeqUL;
     /*One base of last base*/

   if(refLastAlnBaseL > 0)
      alnST->refEndAlnUL = refLastAlnBaseL;
   else alnST->refEndAlnUL = refST->lenSeqUL;
     /*One base of last base*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Add in the query alignment
   \******************************************************/

   #ifdef HIRSCHTWOBIT
      twoBitMvXElmFromStart(qryAlignment, 0);
   #endif
      
   /*The qryAlignment array starts at qryST->offset*/
   qryIndexUL = qryST->offsetUL;
   qryAlnStr += qryIndexUL;

   #ifdef HIRSCHTWOBIT
      bitUC = getTwoBitElm(qryAlignment);
   #else
      bitUC = (uint8_t) *qryAlignment;
   #endif

   while(bitUC != defEndAlnFlag)
   { /*Loop: Add query aligned bases to alnStruct*/
      switch(bitUC)
      { /*Switch: Check the error type*/
         case 0: break;

         case defSnpFlag:
         /*Case: snp*/
            if(qryFirstAlnBaseL < 0) /*Start of alignment*/
               qryFirstAlnBaseL = qryIndexUL;

            qryLastAlnBaseL=qryIndexUL;/*end of alignment*/
            break;
         /*Case: snp*/

         case defMatchFlag:
         /*Case: matches*/
            if(qryFirstAlnBaseL < 0) /*Start of alignment*/
               qryFirstAlnBaseL = qryIndexUL;

            qryLastAlnBaseL=qryIndexUL;/*end of alignment*/
            break;
         /*Case: matches*/

         case defGapFlag:
         /*Case: Insertions*/
            ++numInssL;
            break;
         /*Case: inerstions*/
      } /*Switch: Check the error type*/

      *qryAlnStr = bitUC;

      #ifdef HIRSCHTWOBIT
         twoBitMvToNextElm(qryAlignment);
         bitUC = getTwoBitElm(qryAlignment);
      #else
         ++qryAlignment;
         bitUC = (uint8_t) *qryAlignment;
      #endif

      ++qryIndexUL;
      ++qryAlnStr;
   } /*Loop: Add query aligned bases to alnStruct*/

   if(qryFirstAlnBaseL >= 0) /*Start of alignment*/
      alnST->qryStartAlnUL = qryFirstAlnBaseL;
   else alnST->qryStartAlnUL = qryST->lenSeqUL;
     /*One base of last base*/

   if(qryLastAlnBaseL > 0)
      alnST->qryEndAlnUL = qryLastAlnBaseL;
   else alnST->qryEndAlnUL = qryST->lenSeqUL;
     /*One base of last base*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Add in the alignment stats
   \******************************************************/

   alnST->numInssUL = (unsigned long) numInssL;
   alnST->numDelsUL = (unsigned long) numDelsL;
   alnST->numSnpsUL = (unsigned long) numSnpsL;
   alnST->numMatchesUL = (unsigned long) numMatchesL;

   alnST->lenAlnUL =
      (unsigned long)
      (numInssL + numDelsL + numSnpsL + numMatchesL);

   alnST->refLenUL = refST->lenSeqUL;
   alnST->qryLenUL = qryST->lenSeqUL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun-06 Sec-04:
   ^  - Add in softmasking
   ^  o fun-06 sec-03 sub-01:
   ^    - Add soft masking to the un-aligned ending bases
   ^  o fun-06 sec-03 sub-02:
   ^    - Add soft masking to the un-aligned starting bases
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-01:
   *  - Add soft masking to the un-aligned ending bases
   \******************************************************/

    while(refIndexUL < refST->lenSeqUL)
    { /*Loop: Apply mask to starting reference bases*/
       *refAlnStr = defSoftMaskFlag;
       ++refAlnStr;
       ++refIndexUL;
    } /*Loop: Apply mask to starting reference bases*/

    while(qryIndexUL < qryST->lenSeqUL)
    { /*Loop: Apply mask to starting query bases*/
       *qryAlnStr = defSoftMaskFlag;
       ++qryAlnStr;
       ++qryIndexUL;
    } /*Loop: Apply mask to starting query bases*/

   /******************************************************\
   * Fun-06 Sec-03 Sub-02:
   *  - Add soft masking to the un-aligned starting bases
   \******************************************************/

    refAlnStr = alnST->refAlnStr; 
    while(*refAlnStr == 0)
    { /*Loop: Apply mask to starting reference bases*/
       *refAlnStr = defSoftMaskFlag;
       ++refAlnStr;
    } /*Loop: Apply mask to starting reference bases*/

    qryAlnStr = alnST->qryAlnStr; 
    while(*qryAlnStr == 0)
    { /*Loop: Apply mask to starting query bases*/
       *qryAlnStr = defSoftMaskFlag;
       ++qryAlnStr;
    } /*Loop: Apply mask to starting query bases*/

    return alnST;
} /*twoBitAlnToAlnST*/

/*--------------------------------------------------------\
| Name: vectAddGaps
| Use:
|  - Adds gaps into a array of characters. Uses vectors
|    when large enough.
| Input:
|  - alnAryC:
|    o Array to add gaps to
|  - lenGapL:
|    o number of gaps to add
|  - gapC:
|    o Symbol of gap to add
| Output:
|  - Modifies:
|    o alnAryC to have lenGapL gaps (gapC)
\--------------------------------------------------------*/
void vectAddGaps(
   char *alnAryC, /*Has alignment array to add gaps to*/
   unsigned long lenGapL, /*Size of gap*/
   char gapC,     /*Symbol for a gap*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun-07 TOC: vectAddGaps
   '  - Adds gaps to an array (by vectors if large enough)
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

     long lGap = 0;
     vectI8 gapVectS8;

     if(lenGapL >= defNum8BitElms)
     { /*If: I can add gaps by vectors*/
        gapVectS8 = set1_I8_retVectI8(gapC);

        while(lGap + (defVectBytes - 1) < lenGapL)
        { /*Loop: add in gaps*/
           store_vectI8_retAryI8(alnAryC, gapVectS8);
           alnAryC += defVectBytes;
           lGap += defVectBytes;
        } /*Loop: add in gaps*/
        /*-1 to account for defVectBytes being index 1*/
     } /*If: I can add gaps by vectors*/

     while(lGap < lenGapL)
     { /*Loop: fill in the last gaps*/
       *alnAryC = gapC;
       ++alnAryC;
       ++lGap;
     } /*Loop: fill in the last gaps*/

     return; /*Nothing else to do*/
} /*vectAddGaps*/
