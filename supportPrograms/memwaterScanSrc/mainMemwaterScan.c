/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' mainMemwaterScan SOF: Start Of File
'   - uses a memwater alignment to print the best
'     alignment for each query and reference base
'   o header:
'     - included libraries
'   o fun01: pversion_mainMemwaterScan
'     - prints version number for mainMemwaterScan
'   o fun02: phelp_mainMemwaterScan
'     - prints help message for mainMemwaterScan
'   o fun03: input_mainMemwaterScan
'     - gets user input
'   o main:
'     - driver function to do a waterman alignment
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../genLib/ulCp.h"
#include "../genLib/base10str.h"

#include "../genBio/seqST.h"

#include "../genAln/alnSet.h"
#include "../genAln/indexToCoord.h"
#include "../genAln/memwaterScan.h"

/*.h files only*/
#include "../genLib/endLine.h"
#include "../genAln/alnDefs.h"
#include "../bioTools.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden files
!   o .c  #include "../genLib/fileFun.h"
!   o .c  #include "../genAln/indexToCoord.h"
!   o .h  #include "../genLib/genMath.h" (only .h macros)
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_fqFile_mainMemwaterScan 0
#define def_faFile_mainMemwaterScan 1
#define def_filter_mainMemwaterScan 1/*filter alignments*/

#define def_alnForDir_mainMemwaterScan 1
#define def_alnRevDir_mainMemwaterScan 2
#define def_alnDir_mainMemwaterScan (def_alnForDir_mainMemwaterScan | def_alnRevDir_mainMemwaterScan)

#define def_minPercScore_mainMemwaterScan 0.90f
#define def_minScore_mainMemwaterScan 0

/*-------------------------------------------------------\
| Fun01: pversion_mainMemwaterScan
|   - prints version number for mainMemwaterScan
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints:
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_mainMemwaterScan(
   void *outFILE
){
   fprintf(
     (FILE *) outFILE,
     "memwaterScan from bioTools version: %i-%02i-%02i%s",
     def_year_bioTools,
     def_month_bioTools,
     def_day_bioTools,
     str_endLine
   );
} /*pversion_mainMemwaterScan*/

/*-------------------------------------------------------\
| Fun02: phelp_mainMemwaterScan
|   - prints help message for mainMemwaterScan
| Input:
|   - outFILE:
|     o file to print help message to
| Output:
|   - Prints:
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_mainMemwaterScan(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints help message for mainMemwaterScan
   '   o fun02 sec01:
   '     - print usage block
   '   o fun02 sec02:
   '     - print input block
   '   o fun02 sec03:
   '     - print output block
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec01:
   ^   - print usage block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
     (FILE *) outFILE,
     "mainMemwaterScan -qry query.fa -ref ref.fa"
   );
   fprintf((FILE *) outFILE, " > out.aln%s", str_endLine);

   fprintf(
      (FILE *) outFILE,
      "  - finds starts and ends of possible alignments%s",
      str_endLine
   );

   fprintf(
     (FILE *) outFILE,
     "     between query and a reference sequence%s",
     str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  - WARNING; this system will print out many%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    duplicate results%s",
      str_endLine
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print input block
   ^   o fun02 sec02 sub01:
   ^     - input block header
   ^   o fun02 sec02 sub02:
   ^     - reference sequence options
   ^   o fun02 sec02 sub03:
   ^     - query sequence options
   ^   o fun02 sec02 sub04:
   ^     - output file
   ^   o fun02 sec02 sub05:
   ^     - alignment settings
   ^   o fun02 sec02 sub06:
   ^     - help/version number
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec02 Sub01:
   *   - input block header
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "Input:%s",
      str_endLine
   );

   /*****************************************************\
   * Fun02 Sec02 Sub02:
   *   - reference sequence options
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  -ref reference.fasta: [Required]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o fasta file with reference sequence%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  -ref-fq reference.fastq: [Replaces -ref]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o fastq file with reference sequence%s",
      str_endLine
   );

   /*****************************************************\
   * Fun02 Sec02 Sub03:
   *   - query sequence options
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  -qry reference.fasta: [Required]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o fasta file with query sequence(s)%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  -qry-fq query.fastq: [Replaces -qry]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o fastq file with query sequence(s)%s",
      str_endLine
   );

   /*****************************************************\
   * Fun02 Sec02 Sub04:
   *   - output file
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  -out out.aln: [stdout]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o file to print alignment to%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o use \"-out -\" for stdout%s",
      str_endLine
   );

   /*****************************************************\
   * Fun02 Sec02 Sub05:
   *   - alignment settings
   *   o fun02 sec02 sub05 cat01:
   *     - gap penalties
   *   o fun02 sec02 sub05 cat02:
   *     - score and match matrixes
   *   o fun02 sec02 sub05 cat04:
   *     - minimum score
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat01:
   +   - gap penalties
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      (FILE *) outFILE,
      "  -gap %i: [Optional]%s",
      def_gapOpen_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o gap opening penalty%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  -extend %i: [Optional]%s",
      def_gapExtend_alnDefs,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o gap extension penalty%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat02:
   +   - score and match matrixes
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      (FILE *) outFILE,
      "  -score scoring-matrix.txt: [Optional]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o scoring matrix (alnDefs.h has defaults)%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  -score match-matrix.txt: [Optional]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o match matrix (alnDefs.h has defaults)%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat03:
   +   - minimum percent score
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      (FILE *) outFILE,
      "  -perc-score %.2f: [Optional]%s",
      def_minPercScore_mainMemwaterScan,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum percent of maxium score needed%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "      to print (keep) an alignment%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat04:
   +   - minimum score
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      (FILE *) outFILE,
      "  -score %i: [Optional]%s",
      def_minScore_mainMemwaterScan,
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o minimum score needed to keep alignment%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat05:
   +   - filter alignments
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_filter_mainMemwaterScan)
      fprintf(
         (FILE *) outFILE,
         "  -filter: [yes]%s",
         str_endLine
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -filter: [no]%s",
         str_endLine
      );

   fprintf(
     (FILE *) outFILE,
     "    o filter alignments to remove nested alignments"
   );
   fprintf((FILE *) outFILE, "%s", str_endLine);

   fprintf(
      (FILE *) outFILE,
      "    o disable with `-no-filter`%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat06:
   +   - alignment direction
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(
        def_alnDir_mainMemwaterScan
      & def_alnForDir_mainMemwaterScan
   ) fprintf(
        (FILE *) outFILE,
        "  -forward: [yes]%s",
        str_endLine
     );

   else
      fprintf(
         (FILE *) outFILE,
         "  -forward: [no]%s",
         str_endLine
      );
   fprintf(
     (FILE *) outFILE,
     "    o align query in forward direction%s",
     str_endLine
   );
   fprintf(
     (FILE *) outFILE,
     "    o disable with `-no-forward`%s",
     str_endLine
   );


   if(
        def_alnDir_mainMemwaterScan
      & def_alnRevDir_mainMemwaterScan
   ) fprintf(
        (FILE *) outFILE,
        "  -reverse: [yes]%s",
        str_endLine
     );

   else
      fprintf(
         (FILE *) outFILE,
         "  -reverse: [no]%s",
         str_endLine
      );
   fprintf(
     (FILE *) outFILE,
     "    o align query in reverse direction%s",
     str_endLine
   );
   fprintf(
     (FILE *) outFILE,
     "    o disable with `-no-reverse`%s",
     str_endLine
   );

   /*****************************************************\
   * Fun02 Sec02 Sub06:
   *   - help/version number
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "  -h: print this help message and exit%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  -v: print version number and exit%s",
      str_endLine
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec03:
   ^   - print output block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Output:%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  - prints alignment(s) to -out file or stdout%s",
      str_endLine
   );
} /*phelp_mainMemwaterScan*/

/*-------------------------------------------------------\
| Fun03: input_mainMemwaterScan
|   - gets user input
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o array of c-strings with user input
|   - refFileStrPtr:
|     o c-string pionter to point to reference file name
|   - refTypeSCPtr:
|     o pointer to signed char to hold if reference was
|       a fasta or fastq file
|   - qryFileStrPtr:
|     o c-string pionter to point to query file name
|   - qryTypeSCPtr:
|     o pointer to signed char to hold if query was
|       a fasta or fastq file
|   - outFileStrPtr:
|     o c-string pionter to point to output file name
|   - minPercScoreFPtr:
|     o float pointer to get the minimum percent score
|   - minScoreSLPtr:
|     o signed long pointer to get the minimum score
|   - filterBlPtr:
|     o signed char pointer to get if filtering (1) or
|       not (0)
|   - dirSCPtr:
|     o signed char pointer to get direction of the
|       alignment
|       * def_alnForDir_memwaterScan for foward
|       * def_alnRevDir_memwaterScan for reverse
|       * def_alnForDir_memwaterScan |
|         def_alnRevDir_memwaterScan for both
|   - alnSetSTPtr:
|     o pointer to alnSet struct with alingment settings
| Output:
|   - Modifies:
|     o all input variables to hold/point to user input
|   - Returns:
|     o 0 for no errors
|     o 1 for help message/version number
|     o 2 for unkown input
\-------------------------------------------------------*/
signed char
input_mainMemwaterScan(
   int numArgsSI,
   char *argAryStr[],
   signed char **refFileStrPtr,
   signed char *refTypeSCPtr,
   signed char **qryFileStrPtr,
   signed char *qryTypeSCPtr,
   signed char **outFileStrPtr,
   float *minPercScoreFPtr,
   signed long *minScoreSLPtr,
   signed char *filterBlPtr,
   signed char *dirSCPtr,
   struct alnSet *alnSetSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - gets user input
   '   o fun03 sec01:
   '     - variable declerations
   '   o fun03 sec02:
   '     - check if user input something
   '   o fun03 sec03:
   '     - get user input
   '   o fun03 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned long errUL = 0;
   signed char errSC = 0;
   signed char *errStr = 0;
   signed int siArg = 1;

   FILE *matrixFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - check if user input something
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numArgsSI == 1)
   { /*If: nothing input*/
      phelp_mainMemwaterScan(stdout);
      goto phelp_fun03_sec04;
   } /*If: nothing input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - get user input
   ^   o fun03 sec03 sub01:
   ^     - start loop and check file input
   ^   o fun03 sec03 sub02:
   ^     - check alignment settings
   ^   o fun03 sec03 sub03:
   ^     - check if help message requested
   ^   o fun03 sec03 sub04:
   ^     - check if version number print
   ^   o fun03 sec03 sub05:
   ^     - invalid entry
   ^   o fun03 sec03 sub06:
   ^     - move to next entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec03 Sub01:
   *   - start loop and check file input
   \*****************************************************/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/
      if(
         ! eqlNull_ulCp(
            (signed char *) "-ref",
            (signed char *) argAryStr[siArg]
         )
      ){ /*If: is reference (as fasta)*/
         ++siArg;
         *refFileStrPtr =
            (signed char *) argAryStr[siArg];
         *refTypeSCPtr = def_faFile_mainMemwaterScan;
      } /*If: is reference (as fasta)*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-ref-fq",
            (signed char *) argAryStr[siArg]
         )
      ){ /*If: is reference (as fastq)*/
         ++siArg;
         *refFileStrPtr =
            (signed char *) argAryStr[siArg];
         *refTypeSCPtr = def_fqFile_mainMemwaterScan;
      } /*If: is reference (as fastq)*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-qry",
            (signed char *) argAryStr[siArg]
         )
      ){ /*If: is query (as fasta)*/
         ++siArg;
         *qryFileStrPtr =
            (signed char *) argAryStr[siArg];
         *qryTypeSCPtr = def_faFile_mainMemwaterScan;
      } /*If: is query (as fasta)*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-qry-fq",
            (signed char *) argAryStr[siArg]
         )
      ){ /*If: is query (as fastq)*/
         ++siArg;
         *qryFileStrPtr =
            (signed char *) argAryStr[siArg];
         *qryTypeSCPtr = def_fqFile_mainMemwaterScan;
      } /*If: is query (as fastq)*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-out",
            (signed char *) argAryStr[siArg]
         )
      ){ /*If: is output file*/
         ++siArg;
         *outFileStrPtr =
            (signed char *) argAryStr[siArg];
      } /*If: is output file*/

      /**************************************************\
      * Fun03 Sec03 Sub02:
      *   - check alignment settings
      *   o fun03 sec03 sub02 cat01:
      *     - check general settings
      *   o fun03 sec03 sub02 cat02:
      *     - check scoring matrix input
      *   o fun03 sec03 sub02 cat03:
      *     - check scoring match input
      *   o fun03 sec03 sub02 cat04:
      *     - check minimum percent score
      *   o fun03 sec03 sub02 cat05:
      *     - check minimum score
      *   o fun03 sec03 sub02 cat06:
      *     - check if filtering
      *   o fun03 sec03 sub02 cat07:
      *     - check direction
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat01:
      +   - check general settings
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-gap",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: is gap opening penalty*/
         ++siArg;

         errStr +=
            strToSS_base10str(
               (signed char *) argAryStr[siArg],
               &alnSetSTPtr->gapSS
            );

         if(*errStr)
         { /*If: non-numeric*/
            fprintf(
               stderr,
               "-gap %s is non-numeric%s",
               argAryStr[siArg],
               str_endLine
            );
               
            goto err_fun03_sec04;     
         } /*If: non-numeric*/
      } /*Else If: is gap opening penalty*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-extend",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: is gap extension penalty*/
         ++siArg;

         errStr +=
            strToSS_base10str(
               (signed char *) argAryStr[siArg],
               &alnSetSTPtr->extendSS
            );

         if(*errStr)
         { /*If: non-numeric*/
            fprintf(
               stderr,
               "-extend %s is non-numeric%s",
               argAryStr[siArg],
               str_endLine
            );
               
            goto err_fun03_sec04;     
         } /*If: non-numeric*/
      } /*Else If: is gap extension penalty*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat02:
      +   - check scoring matrix input
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-score",
            (signed char *) argAryStr[siArg]
         )
      ){ /*If: is scoring matrix*/
         ++siArg;

         matrixFILE =
            fopen(
               argAryStr[siArg],
               "r"
           );

         if(! matrixFILE)
         { /*If: could not open scoring matrix*/
            fprintf(
               stderr,
               "unable to open -score %s%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*If: could not open scoring matrix*/

         errUL =
            readScoreFile_alnSet(
               alnSetSTPtr,
               matrixFILE
            );

         fclose(matrixFILE);
         matrixFILE = 0;

         if(errUL)
         { /*If: invalid entry*/
            fprintf(
               stderr,
               "entry %lu in -score %s is invalid%s",
               errUL,
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*If: invalid entry*/
      } /*If: is scoring matrix*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat03:
      +   - check match matrix input
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-match",
            (signed char *) argAryStr[siArg]
         )
      ){ /*If: is match matrix*/
         ++siArg;

         matrixFILE =
            fopen(
               argAryStr[siArg],
               "r"
           );

         if(! matrixFILE)
         { /*If: could not open scoring matrix*/
            fprintf(
               stderr,
               "unable to open -match %s%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*If: could not open scoring matrix*/

         errUL =
            readMatchFile_alnSet(
               alnSetSTPtr,
               matrixFILE
            );

         fclose(matrixFILE);
         matrixFILE = 0;

         if(errUL)
         { /*If: invalid entry*/
            fprintf(
               stderr,
               "entry %lu in -match %s is invalid%s",
               errUL,
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*If: invalid entry*/
      } /*If: is match matrix*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat04:
      +   - check minimum percent score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-perc-score",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: minimum percent score*/
         ++siArg;
         errStr = (signed char *) argAryStr[siArg];
         errStr +=
            strToF_base10str(errStr, minPercScoreFPtr);

         if(*errStr > 32)
         { /*If: invalid input*/
            fprintf(
               stderr,
               "-perc-score %s is non-numeric%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*If: invalid input*/

         else if(
               *minPercScoreFPtr > 1
            && *minPercScoreFPtr <= 100
         ) *minPercScoreFPtr /= 100; /*percent input*/

         else if(*minPercScoreFPtr > 100)
         { /*Else If: value greater then 100%*/
            fprintf(
               stderr,
               "-perc-score %s can not be over 100%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*Else If: value greater then 100%*/

         else if(*minPercScoreFPtr < 0)
         { /*Else If: negative input*/
            fprintf(
               stderr,
               "-perc-score %s can not be negative%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*Else If: negative input*/
      }  /*Else If: minimum percent score*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat05:
      +   - check minimum score
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-score",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: minimum score*/
         ++siArg;
         errStr = (signed char *) argAryStr[siArg];
         errStr +=
            strToSL_base10str(errStr, minScoreSLPtr);

         if(*errStr > 32)
         { /*If: invalid input*/
            fprintf(
               stderr,
               "-score %s is non-numeric%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*If: invalid input*/

         else if(*minPercScoreFPtr < 0)
         { /*Else If: negative input*/
            fprintf(
               stderr,
               "-score %s can not be negative%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec04;     
         } /*Else If: negative input*/
      }  /*Else If: minimum score*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat06:
      +   - check if filtering
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-filter",
            (signed char *) argAryStr[siArg]
         )
      ) *filterBlPtr = 1;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-filter",
            (signed char *) argAryStr[siArg]
         )
      ) *filterBlPtr = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat07:
      +   - check direction
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-forward",
            (signed char *) argAryStr[siArg]
         )
      ) *dirSCPtr |= def_alnForDir_mainMemwaterScan;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-forward",
            (signed char *) argAryStr[siArg]
         )
      ) *dirSCPtr &= ~def_alnForDir_mainMemwaterScan;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-reverse",
            (signed char *) argAryStr[siArg]
         )
      ) *dirSCPtr |= def_alnRevDir_mainMemwaterScan;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-reverse",
            (signed char *) argAryStr[siArg]
         )
      ) *dirSCPtr &= ~def_alnRevDir_mainMemwaterScan;

      /**************************************************\
      * Fun03 Sec03 Sub03:
      *   - check if help message requested
      \**************************************************/
      
      else if(
         ! eqlNull_ulCp(
            (signed char *) "-h",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "--h",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/


      else if(
         ! eqlNull_ulCp(
            (signed char *) "help",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-help",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "--help",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/

      /**************************************************\
      * Fun03 Sec03 Sub04:
      *   - check if version number print
      \**************************************************/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-v",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "--v",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "version",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-version",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "--version",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwaterScan(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      /**************************************************\
      * Fun03 Sec03 Sub05:
      *   - invalid entry
      \**************************************************/

      else
      { /*Else: invalid input*/
         fprintf(
            stderr,
            "%s is not recognized%s",
            argAryStr[siArg],
            str_endLine
         );

         goto err_fun03_sec04;     
      } /*Else: invalid input*/

      /**************************************************\
      * Fun03 Sec03 Sub06:
      *   - move to next entry
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*no errors*/
   errSC = 0;
   goto cleanUp_fun03_sec04;

   phelp_fun03_sec04:;
   errSC = 1;
   goto cleanUp_fun03_sec04;

   pversion_fun03_sec04:;
   errSC = 1;
   goto cleanUp_fun03_sec04;

   err_fun03_sec04:;
   errSC = 2;
   goto cleanUp_fun03_sec04;

   cleanUp_fun03_sec04:;

   if(
         matrixFILE
      && matrixFILE != stdin
      && matrixFILE != stdout
   ) fclose(matrixFILE);

   matrixFILE = 0;
   return errSC;
} /*input_mainMemwaterScan*/

/*-------------------------------------------------------\
| Main:
|   - driver function to do a waterman alignment
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o array of c-strings with user input
| Output:
|   - prints alignment as sam file
\-------------------------------------------------------*/
int
main(
   int numArgsSI,
   char *argAryStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '   - driver function to do a waterman alignment
   '   o main sec01:
   '     - variable declerations
   '   o main sec02:
   '     - initialize structs, get input, and check input
   '   o main sec03:
   '     - align sequences
   '   o main sec04:
   '     - clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   signed long seqSL = 0;   /*query sequence on*/
   signed char filterBl = def_filter_mainMemwaterScan;
   signed char dirFlagSC = def_alnDir_mainMemwaterScan;
   float minPercScoreF=def_minPercScore_mainMemwaterScan;
   signed long minScoreSL = def_minScore_mainMemwaterScan;

   signed char didRevBl = 0;
   signed char *qryFileStr = 0;
   signed char qryTypeSC = def_fqFile_mainMemwaterScan;

   signed char *refFileStr = 0;
   signed char refTypeSC = def_fqFile_mainMemwaterScan;


   /*values for final return*/
   signed char dirCharSC = 0;
   signed char *outFileStr = 0;
   signed int siNt = 0; 

   /*values for final return*/
   signed long maxScoreSL = 0;
   signed long scoreSL = 0;
   float percScoreF = 0;
   float scoreF = 0;

   unsigned long refStartUL = 0;
   unsigned long refEndUL = 0;
   unsigned long qryStartUL = 0;
   unsigned long qryEndUL = 0;


   signed char *tmpStr = 0;
   struct seqST qryStackST;
   struct seqST refStackST;
   struct alnSet setStackST;
   struct aln_memwaterScan alnStackST;

   FILE *seqFILE = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - initialize structures, get input, and check input
   ^   o main sec02 sub01:
   ^     - initialize structures
   ^   o main sec02 sub02:
   ^     - get input
   ^   o main sec02 sub03:
   ^     - get reference sequence
   ^   o main sec02 sub04:
   ^     - get first query sequence
   ^   o main sec02 sub05:
   ^     - open output file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize structures
   \*****************************************************/

   init_seqST(&qryStackST);
   init_seqST(&refStackST);
   init_alnSet(&setStackST);
   init_aln_memwaterScan(&alnStackST);

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get input
   \*****************************************************/

   errSC =
      input_mainMemwaterScan(
         numArgsSI,
         argAryStr,
         &refFileStr,
         &refTypeSC,
         &qryFileStr,
         &qryTypeSC,
         &outFileStr,
         &minPercScoreF,
         &minScoreSL,
         &filterBl,
         &dirFlagSC,
         &setStackST
      );

   if(errSC)
   { /*If: had input error*/
      --errSC;
      goto cleanUp_main_sec04;
   } /*If: had input error*/

   if(! dirFlagSC)
   { /*If: no alignment*/
      fprintf(
         stderr,
         "both -no-forward and -no-reverse are used,%s",
         str_endLine
      );

      fprintf(
         stderr,
         " which prevents any aligmnet, pick one%s",
         str_endLine
      );

      errSC = 8;
      goto cleanUp_main_sec04;
   } /*If: no alignment*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - get reference sequence
   *   o main sec02 sub03 cat01:
   *     - open reference file
   *   o main sec02 sub03 cat02:
   *     - get reference sequence
   *   o main sec02 sub03 cat03:
   *     - deal with errors
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat01:
   +   - open reference file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(! refFileStr)
   { /*If: no reference file input*/
      fprintf(
         stderr,
         "no reference file provied with -ref%s",
         str_endLine
      );

      goto fileErr_main_sec04;
   } /*If: no reference file input*/

   seqFILE = fopen((char *) refFileStr, "r");

   if(! seqFILE)
   { /*If: file error*/
      if(refTypeSC == def_fqFile_mainMemwaterScan)
         fprintf( stderr,
            "could not open -ref-fq %s%s",
            refFileStr,
            str_endLine
         );

      else
         fprintf( stderr,
            "could not open -ref %s%s",
            refFileStr,
            str_endLine
         );

      goto fileErr_main_sec04;
   } /*If: file error*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat02:
   +   - get reference sequence
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(refTypeSC == def_fqFile_mainMemwaterScan)
      errSC = getFq_seqST(seqFILE, &refStackST);
   else
      errSC = getFa_seqST(seqFILE, &refStackST);

   if(! errSC)
   { /*If: no errors; prepare id for priting*/
      /*remove white space from reference id*/
      tmpStr = refStackST.idStr;
      while(*tmpStr++ > 32)
         ;
      --tmpStr;
      *tmpStr = '\0';
   } /*If: no errors; prepare id for priting*/

   fclose(seqFILE);
   seqFILE = 0;

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat03:
   +   - deal with errors
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(
         errSC
      && errSC != def_EOF_seqST
   ){ /*If: had an error*/
      if(errSC & def_memErr_seqST)
         fprintf( stderr, "MEMORY ERROR reading");
      else
         fprintf( stderr, "file error with");

      if(refTypeSC == def_fqFile_mainMemwaterScan)
         fprintf(
           stderr,
           " -ref-fq %s%s",
           refFileStr,
           str_endLine
         );

      else
         fprintf(
           stderr,
           " -ref %s%s",
           refFileStr,
           str_endLine
         );

      if(errSC & def_memErr_seqST)
         goto memErr_main_sec04;
      else
         goto fileErr_main_sec04;
   } /*If: had an error*/

   /*****************************************************\
   * Main Sec02 Sub04:
   *   - get frist query sequence
   *   o main sec02 sub04 cat01:
   *     - open query file
   *   o main sec02 sub04 cat02:
   *     - get query sequence
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub04 Cat01:
   +   - open query file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(! qryFileStr)
   { /*If: no query file input*/
      fprintf(
         stderr,
         "no query file provied with -qry%s",
         str_endLine
      );

      goto fileErr_main_sec04;
   } /*If: no query file input*/

   seqFILE = fopen((char *) qryFileStr, "r");

   if(! seqFILE)
   { /*If: file error*/
      if(qryTypeSC == def_fqFile_mainMemwaterScan)
         fprintf( stderr,
            "could not open -qry-fq %s%s",
            qryFileStr,
            str_endLine
         );

      else
         fprintf( stderr,
            "could not open -qry %s%s",
            qryFileStr,
            str_endLine
         );

      goto fileErr_main_sec04;
   } /*If: file error*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub04 Cat02:
   +   - get query sequence
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(qryTypeSC == def_fqFile_mainMemwaterScan)
      errSC = getFq_seqST(seqFILE, &qryStackST);
   else
      errSC = getFa_seqST(seqFILE, &qryStackST);

   if(
         errSC
      && errSC != def_EOF_seqST
   ){ /*If: had an error*/
      if(errSC & def_memErr_seqST)
         fprintf(stderr, "MEMORY ERROR reading");
      else
         fprintf(stderr, "file error with");

      if(qryTypeSC == def_fqFile_mainMemwaterScan)
         fprintf(stderr, " -qry-fq %s", qryFileStr);
      else
         fprintf(stderr, " -qry %s", qryFileStr);
      fprintf(stderr, "%s", str_endLine);

      if(errSC & def_memErr_seqST)
         goto memErr_main_sec04;
      else
         goto fileErr_main_sec04;
   } /*If: had an error*/

   /*****************************************************\
   * Main Sec02 Sub05:
   *   - open output file
   \*****************************************************/

   if(! outFileStr || *outFileStr == '-')
      outFILE = stdout;

   else
   { /*Else: output file input*/
      outFILE = fopen((char *) outFileStr, "w");

      if(! outFILE)
      { /*If: could not open output file*/
         fprintf(
            stderr,
            "unable to open -out %s%s",
            outFileStr,
            str_endLine
         );

         goto fileErr_main_sec04;
      } /*If: could not open output file*/
   } /*Else: output file input*/
  
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - align sequences
   ^   o main sec03 sub01:
   ^     - setup/start loop and prepare for alignment
   ^   o main sec03 sub02:
   ^     - align and print out alignment
   ^   o main sec03 sub03:
   ^     - get next query sequence
   ^   o main sec03 sub04:
   ^     - finished alingments, check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub01:
   *   - setup/start loop and prepare for alignment
   \*****************************************************/

   /*print out the header*/
   fprintf(
      outFILE,
      "qry_id\tref_id\tdir\tscore\tperc_score\tmax_score"
   );

   fprintf(
      outFILE,
      "\tref_start\tref_end\tref_len\tqry_start\tqry_end"
   );

   fprintf(
      outFILE,
      "\tqry_len%s",
      str_endLine
   );

   seqToIndex_alnSet(refStackST.seqStr);

   do
   { /*Loop: align all query sequences*/
      ++seqSL;

      /*remove white space from query id*/
      tmpStr = qryStackST.idStr;
      while(*tmpStr++ > 32) ;
      --tmpStr;
      *tmpStr = '\0';

      seqToIndex_alnSet(qryStackST.seqStr);

      if(! (dirFlagSC & def_alnForDir_mainMemwaterScan) )
         goto revAln_main_sec03_sub03;
      else
         dirCharSC = 'F';

      algin_main_sec03_sub01:;
      maxScoreSL =
         maxScore_alnSet(
            qryStackST.seqStr,
            0,
            qryStackST.seqLenSL,
            &setStackST
         );

      /**************************************************\
      * Main Sec03 Sub02:
      *   - align and print out alignment
      *   o main sec03 sub02 cat01:
      *     - align input
      *   o main sec03 sub02 cat02:
      *     - print the best score for each base
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec03 Sub02 Cat01:
      +   - align input
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      scoreSL =
         memwaterScan(
            &qryStackST,
            &refStackST,
            &alnStackST,
            &setStackST
         );

      if(scoreSL < 0)
      { /*If: memory error*/
         errSC = def_memErr_seqST;
         break;
      } /*If: memory error*/

      if(filterBl)
      { /*If: filtering alignments*/
         filter_memwaterScan(
            &alnStackST,
            maxScoreSL * minPercScoreF, /*percent cutoff*/
            minScoreSL                  /*hard cutoff*/
         );
      } /*If: filtering alignments*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Main Sec03 Sub02 Cat02:
      +   - print the best score for each base
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      for(siNt = 0; siNt < alnStackST.outLenSL; ++siNt)
      { /*Loop: print out best scores*/
         if(alnStackST.scoreArySL[siNt] < minScoreSL)
            continue; /*under minimum score*/
 
         scoreF = alnStackST.scoreArySL[siNt];
         if(scoreF < 0)
            continue; /*no score assigned, very unlikely*/
         scoreF /= (float) def_scoreAdj_alnDefs;

         percScoreF =
              (float) alnStackST.scoreArySL[siNt]
            / (float) maxScoreSL;
         if(percScoreF < minPercScoreF)
            continue; /*score is to low*/

         indexToCoord(
            alnStackST.refLenSI,
            alnStackST.startArySL[siNt],
            &refStartUL,
            &qryStartUL
         );
         indexToCoord(
            alnStackST.refLenSI,
            alnStackST.endArySL[siNt],
            &refEndUL,
            &qryEndUL
         );

         fprintf(
           outFILE,
           "%s\t%s\t%c\t%0.2f\t%0.2f\t%li",
           qryStackST.idStr,
           refStackST.idStr,
           dirCharSC,
           scoreF,
           percScoreF,
           maxScoreSL
         ); /*non-coodinates*/
         fprintf(
            outFILE,
            "\t%lu\t%lu\t%li\t%lu\t%lu\t%li%s",
            refStartUL + 1,
            refEndUL + 1,
            refStackST.seqLenSL,
            qryStartUL + 1,
            qryEndUL + 1,
            qryStackST.seqLenSL,
            str_endLine
         );
      } /*Loop: print out best scores*/

      /**************************************************\
      * Main Sec03 Sub03:
      *   - get next query sequence
      \**************************************************/

      revAln_main_sec03_sub03:;
         if(didRevBl)
            didRevBl = 0;
         else if(
            dirFlagSC & def_alnRevDir_mainMemwaterScan
         ){ /*Else If: doing a reverse alignment*/
            revCmpIndex_alnSet(
               qryStackST.seqStr,
               qryStackST.qStr,
               qryStackST.seqLenSL
            );

            dirCharSC = 'R';
            didRevBl = 1;
            goto algin_main_sec03_sub01;
         } /*Else If: doing a reverse alignment*/

      if(qryTypeSC == def_fqFile_mainMemwaterScan)
         errSC = getFq_seqST(seqFILE, &qryStackST);
      else
         errSC = getFa_seqST(seqFILE, &qryStackST);
   }  while(! errSC);
      /*Loop: align all query sequences*/

   /*****************************************************\
   * Main Sec03 Sub04:
   *   - finished alingments, check for errors
   \*****************************************************/

   if(errSC != def_EOF_seqST)
   { /*If: memory error*/
      fprintf(
         stderr,
         "MEMORY ERROR on %li%s",
         seqSL,
         str_endLine
      );

      goto memErr_main_sec04;
   } /*If: memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto cleanUp_main_sec04;


   memErr_main_sec04:;
      errSC = 1;
      goto cleanUp_main_sec04;
   
   fileErr_main_sec04:;
      errSC = 2;
      goto cleanUp_main_sec04;

   cleanUp_main_sec04:;
      freeStack_seqST(&qryStackST);
      freeStack_seqST(&refStackST);
      freeStack_alnSet(&setStackST);
      freeStack_aln_memwaterScan(&alnStackST);
   
      if(! seqFILE) ;
      else if(seqFILE == stdin) ;
      else if(seqFILE == stdout) ;
      else if(seqFILE == stderr) ;
      else fclose(seqFILE);
   
      seqFILE = 0;
   
      if(! outFILE) ;
      else if(outFILE == stdin) ;
      else if(outFILE == stdout) ;
      else if(outFILE == stderr) ;
      else fclose(outFILE);
   
      outFILE = 0;
   
      return errSC;
} /*main*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible, this code is under the
:   MIT license.
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
