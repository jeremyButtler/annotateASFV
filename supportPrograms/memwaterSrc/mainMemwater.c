/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' mainMemwater SOF: Start Of File
'   - does a using a memwater alingment
'   o header:
'     - included libraries
'   o fun01: pversion_mainMemwater
'     - prints version number for mainMemwater
'   o fun02: phelp_mainMemwater
'     - prints help message for mainMemwater
'   o fun03: input_mainMemwater
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

#include "../genLib/charCp.h"
#include "../genLib/base10str.h"

#include "../genBio/seqST.h"
#include "../genBio/kmerCnt.h"

#include "../genAln/alnSet.h"
#include "../genAln/memwater.h"

/*.h files only*/
#include "../genLib/endLine.h"
#include "../genAln/alnDefs.h"
#include "../bioTools.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden files
!   o .c  #include "../genLib/ulCp.h"
!   o .c  #include "../genLib/numToStr.h"
!   o .c  #include "../genLib/shellSort.h"
!   o .c  #include "../genLib/fileFun.h"
!   o .c  #include "../genAln/indexToCoord.h"
!   o .h  #include "../genLib/genMath.h" (only .h macros)
!   o .h  #include "../genBio/ntTo2Bit.h"
!   o .h  #include "../genBio/ntTo5Bit.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_fqFile_mainMemwater 0
#define def_faFile_mainMemwater 1

#define def_lenKmer_mainMemwater 5

#define def_kmerScan_mainMemwater 0
   /*look for direction with kmer counting*/

/*-------------------------------------------------------\
| Fun01: pversion_mainMemwater
|   - prints version number for mainMemwater
| Input:
|   - outFILE:
|     o file to print version number to
| Output:
|   - Prints:
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_mainMemwater(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "memwater from bioTools version: %i-%02i-%02i%s",
      def_year_bioTools,
      def_month_bioTools,
      def_day_bioTools,
      str_endLine
   );
} /*pversion_mainMemwater*/

/*-------------------------------------------------------\
| Fun02: phelp_mainMemwater
|   - prints help message for mainMemwater
| Input:
|   - outFILE:
|     o file to print help message to
| Output:
|   - Prints:
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_mainMemwater(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints help message for mainMemwater
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
     "mainMemwater -qry query.fa -ref ref.fa > out.aln%s",
     str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  - does memwater alignment on two sequences%s",
      str_endLine
   );

   fprintf(
     (FILE *) outFILE,
     "  - only finds score + start and end coordinates%s",
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
   *     - kmer scanning
   *   o fun02 sec02 sub05 cat01:
   *     - gap penalties
   *   o fun02 sec02 sub05 cat02:
   *     - score and match matrixes
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat01:
   +   - kmer scanning
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_kmerScan_mainMemwater)
      fprintf(
         (FILE *) outFILE,
         "  -scan: [Optional; Yes]%s",
         str_endLine
      );

   else
      fprintf(
         (FILE *) outFILE,
         "  -scan: [Optional; No]%s",
         str_endLine
      );

   fprintf(
      (FILE *) outFILE,
      "    o scan for direction using kmers%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o reference and query need to be similar%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "      lengths%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "    o disable with \"-no-scan\"%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun02 Sec02 Sub05 Cat02:
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
   + Fun02 Sec02 Sub05 Cat03:
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
} /*phelp_mainMemwater*/

/*-------------------------------------------------------\
| Fun03: input_mainMemwater
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
|   - kmerScanBl:
|     o 1: use kmer counting to find reference/direction
|     o 0: assume user knows what they are doing
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
input_mainMemwater(
   int numArgsSI,
   char *argAryStr[],
   signed char **refFileStrPtr,
   signed char *refTypeSCPtr,
   signed char **qryFileStrPtr,
   signed char *qryTypeSCPtr,
   signed char **outFileStrPtr,
   signed char *kmerScanBl,
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
      phelp_mainMemwater(stdout);
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
         ! eql_charCp(
            (signed char *) "-ref",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*If: is reference (as fasta)*/
         ++siArg;
         *refFileStrPtr =
            (signed char *) argAryStr[siArg];
         *refTypeSCPtr = def_faFile_mainMemwater;
      } /*If: is reference (as fasta)*/

      else if(
         ! eql_charCp(
            (signed char *) "-ref-fq",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*If: is reference (as fastq)*/
         ++siArg;
         *refFileStrPtr =
            (signed char *) argAryStr[siArg];
         *refTypeSCPtr = def_fqFile_mainMemwater;
      } /*If: is reference (as fastq)*/

      else if(
         ! eql_charCp(
            (signed char *) "-qry",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*If: is query (as fasta)*/
         ++siArg;
         *qryFileStrPtr =
            (signed char *) argAryStr[siArg];
         *qryTypeSCPtr = def_faFile_mainMemwater;
      } /*If: is query (as fasta)*/

      else if(
         ! eql_charCp(
            (signed char *) "-qry-fq",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*If: is query (as fastq)*/
         ++siArg;
         *qryFileStrPtr =
            (signed char *) argAryStr[siArg];
         *qryTypeSCPtr = def_fqFile_mainMemwater;
      } /*If: is query (as fastq)*/

      else if(
         ! eql_charCp(
            (signed char *) "-out",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
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
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec03 Sub02 Cat01:
      +   - check general settings
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eql_charCp(
            (signed char *) "-scan",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ) *kmerScanBl = 1;

      else if(
         ! eql_charCp(
            (signed char *) "-no-scan",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ) *kmerScanBl = 0;

      else if(
         ! eql_charCp(
            (signed char *) "-gap",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
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
         ! eql_charCp(
            (signed char *) "-extend",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
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
         ! eql_charCp(
            (signed char *) "-score",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
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
         ! eql_charCp(
            (signed char *) "-match",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
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

      /**************************************************\
      * Fun03 Sec03 Sub03:
      *   - check if help message requested
      \**************************************************/
      
      else if(
         ! eql_charCp(
            (signed char *) "-h",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eql_charCp(
            (signed char *) "--h",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/


      else if(
         ! eql_charCp(
            (signed char *) "help",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eql_charCp(
            (signed char *) "-help",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/

      else if(
         ! eql_charCp(
            (signed char *) "--help",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: help message requested*/
         phelp_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: help message requested*/

      /**************************************************\
      * Fun03 Sec03 Sub04:
      *   - check if version number print
      \**************************************************/

      else if(
         ! eql_charCp(
            (signed char *) "-v",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eql_charCp(
            (signed char *) "--v",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eql_charCp(
            (signed char *) "version",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eql_charCp(
            (signed char *) "-version",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwater(stdout);
         goto pversion_fun03_sec04;
      } /*Else If: version number requested*/

      else if(
         ! eql_charCp(
            (signed char *) "--version",
            (signed char *) argAryStr[siArg],
            (signed char) '\0'
         )
      ){ /*Else If: version number requested*/
         pversion_mainMemwater(stdout);
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
} /*input_mainMemwater*/

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
   signed long scoreSL = 0;
   signed long seqSL = 0;   /*query sequence on*/

   signed long refStartSL = 0;
   signed long refEndSL = 0;

   signed long qryStartSL = 0;
   signed long qryEndSL = 0;

   signed char *qryFileStr = 0;
   signed char qryTypeSC = def_fqFile_mainMemwater;

   signed char *refFileStr = 0;
   signed char refTypeSC = def_fqFile_mainMemwater;

   signed char kmerScanBl = def_kmerScan_mainMemwater;
     /*scan for direction using kmers*/

   signed char *tmpStr = 0;

   signed char *outFileStr = 0;

   struct seqST qryStackST;
   struct seqST *refSTPtr = 0;
   struct alnSet setStackST;

   struct kmerCnt kmerStackST;
   unsigned char lenKmerUC = def_lenKmer_mainMemwater;
   unsigned int maxKmerUI = 0;
   signed int numKmersSI = 0;
   signed int *kmerHeapArySI = 0;
   signed int *cntHeapArySI = 0;

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
   ^     - allocate memory for kmer arrays
   ^   o main sec02 sub04:
   ^     - get first reference sequence
   ^   o main sec02 sub05:
   ^     - get first query sequence
   ^   o main sec02 sub06:
   ^     - open output file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize structures
   \*****************************************************/

   init_seqST(&qryStackST);
   init_alnSet(&setStackST);
   init_kmerCnt(&kmerStackST);

   maxKmerUI = 1;

   for(
      scoreSL = 0;
      scoreSL < lenKmerUC;
      ++scoreSL
   ) maxKmerUI <<= 2; /*find max number kmers in table*/

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get input
   \*****************************************************/

   errSC =
      input_mainMemwater(
         numArgsSI,
         argAryStr,
         &refFileStr,
         &refTypeSC,
         &qryFileStr,
         &qryTypeSC,
         &outFileStr,
         &kmerScanBl,
         &setStackST
      );

   if(errSC)
   { /*If: had input error*/
      --errSC;
      goto cleanUp_main_sec04_sub04;
   } /*If: had input error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - allocate memory for kmer arrays
   \*****************************************************/

   kmerHeapArySI =
      malloc((maxKmerUI + 1) * sizeof(signed int));

   if(! kmerHeapArySI)
      goto memErr_main_sec04_sub02;

   cntHeapArySI =
      malloc((maxKmerUI + 1) * sizeof(signed int));

   if(! cntHeapArySI)
      goto memErr_main_sec04_sub02;

   /*****************************************************\
   * Main Sec02 Sub04:
   *   - get reference sequence
   *   o main sec02 sub04 cat01:
   *     - open reference file
   *   o main sec02 sub04 cat02:
   *     - get reference sequence
   *   o main sec02 sub04 cat03:
   *     - deal with errors
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub04 Cat01:
   +   - open reference file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(! refFileStr)
   { /*If: no reference file input*/
      fprintf(
         stderr,
         "no reference file provied with -ref%s",
         str_endLine
      );

      goto fileErr_main_sec04_sub03;
   } /*If: no reference file input*/

   seqFILE =
      fopen(
         (char *) refFileStr,
         "r"
      );

   if(! seqFILE)
   { /*If: file error*/
      if(refTypeSC == def_fqFile_mainMemwater)
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

      goto fileErr_main_sec04_sub03;
   } /*If: file error*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub04 Cat02:
   +   - get reference sequence
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   setup_kmerCnt(
      &kmerStackST,
      lenKmerUC
   );

   if(refTypeSC == def_fqFile_mainMemwater)
      errSC = getFq_seqST(seqFILE, kmerStackST.forSeqST);

   else
      errSC = getFa_seqST(seqFILE, kmerStackST.forSeqST);

   if(! errSC)
   { /*If: no errors; prepare kmerCnt structure*/
      /*remove white space from reference id*/
      tmpStr = kmerStackST.forSeqST->idStr;
      while(*tmpStr++ > 32) ;
      --tmpStr;
      *tmpStr = '\0';

      if(kmerScanBl)
         errSC =
            addSeq_kmerCnt(
               &kmerStackST,
               kmerStackST.forSeqST
            );
   } /*If: no errors; prepare kmerCnt structure*/

   fclose(seqFILE);
   seqFILE = 0;

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub04 Cat03:
   +   - deal with errors
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(
         errSC
      && errSC != def_EOF_seqST
   ){ /*If: had an error*/
      if(errSC & def_memErr_seqST)
         fprintf(
            stderr,
            "MEMORY ERROR reading"
         );

      else
         fprintf(
            stderr,
            "file error with"
         );

      if(refTypeSC == def_fqFile_mainMemwater)
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
         goto memErr_main_sec04_sub02;
      else
         goto fileErr_main_sec04_sub03;
   } /*If: had an error*/

   /*****************************************************\
   * Main Sec02 Sub05:
   *   - get frist query sequence
   *   o main sec02 sub03 cat01:
   *     - open query file
   *   o main sec02 sub03 cat02:
   *     - get reference sequence
   *   o main sec02 sub03 cat03:
   *     - deal with errors
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat01:
   +   - open query file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(! qryFileStr)
   { /*If: no query file input*/
      fprintf(
         stderr,
         "no query file provied with -qry%s",
         str_endLine
      );

      goto fileErr_main_sec04_sub03;
   } /*If: no query file input*/

   seqFILE =
      fopen(
         (char *) qryFileStr,
         "r"
      );

   if(! seqFILE)
   { /*If: file error*/
      if(qryTypeSC == def_fqFile_mainMemwater)
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

      goto fileErr_main_sec04_sub03;
   } /*If: file error*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat02:
   +   - get query sequence
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(qryTypeSC == def_fqFile_mainMemwater)
      errSC = getFq_seqST(seqFILE, &qryStackST);

   else
      errSC = getFa_seqST(seqFILE, &qryStackST);

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub03 Cat03:
   +   - deal with errors
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(
         errSC
      && errSC != def_EOF_seqST
   ){ /*If: had an error*/
      if(errSC & def_memErr_seqST)
         fprintf(stderr, "MEMORY ERROR reading");
      else
         fprintf(stderr, "file error with");

      if(qryTypeSC == def_fqFile_mainMemwater)
         fprintf(stderr, " -qry-fq %s", qryFileStr);
      else
         fprintf(stderr, " -qry %s", qryFileStr);
      fprintf(stderr, "%s", str_endLine);

      if(errSC & def_memErr_seqST)
         goto memErr_main_sec04_sub02;
      else
         goto fileErr_main_sec04_sub03;
   } /*If: had an error*/

   /*****************************************************\
   * Main Sec02 Sub06:
   *   - open output file
   \*****************************************************/

   if(! outFileStr || *outFileStr == '-')
      outFILE = stdout;

   else
   { /*Else: output file input*/
      outFILE =
         fopen(
            (char *) outFileStr,
            "w"
         );

      if(! outFILE)
      { /*If: could not open output file*/
         fprintf(
            stderr,
            "unable to open -out %s%s",
            outFileStr,
            str_endLine
         );

         goto fileErr_main_sec04_sub03;
      } /*If: could not open output file*/
   } /*Else: output file input*/
  
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^   - align sequences
   ^   o main sec03 sub01:
   ^     - setup/start loop and prepare for alignment
   ^   o main sec03 sub02:
   ^     - scan for direction (using kmer counting)
   ^   o main sec03 sub03:
   ^     - align and print out alignment
   ^   o main sec03 sub04:
   ^     - get next query sequence
   ^   o main sec03 sub05:
   ^     - finished alingments, check for errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub01:
   *   - setup/start loop and prepare for alignment
   \*****************************************************/

   /*print out the header*/
   fprintf(
      outFILE,
      "qry_id\tref_id\tscore\tdir\tref_start\tref_end"
   );

   fprintf(
      outFILE,
      "\tqry_start\tqry_end%s",
      str_endLine
   );

   seqToIndex_alnSet(kmerStackST.forSeqST->seqStr);

   if(kmerScanBl)
      seqToIndex_alnSet(kmerStackST.revSeqST->seqStr);
   else
       refSTPtr = kmerStackST.forSeqST;

   do
   { /*Loop: align all query sequences*/
      ++seqSL;

      /*remove white space from query id*/
      tmpStr = qryStackST.idStr;
      while(*tmpStr++ > 32) ;
      --tmpStr;
      *tmpStr = '\0';

      /**************************************************\
      * Main Sec03 Sub02:
      *   - scan for direction (using kmer counting)
      \**************************************************/

      if(kmerScanBl)
      { /*If: finding direction with kmers*/
         numKmersSI =
            ntToKmerAry_kmerCnt(
               &qryStackST,
               lenKmerUC,
               kmerHeapArySI,
               cntHeapArySI
            ); /*get kmers in sequence*/

         /*find the reference direction*/
         numKmersSI =
            get_kmerCnt(
               &kmerStackST,
               kmerHeapArySI,
               cntHeapArySI
            ); /*get number kmers shared between ref/qry*/
            
         /*find the direction*/
         if(numKmersSI < 0)
            refSTPtr = kmerStackST.revSeqST;
         else
            refSTPtr = kmerStackST.forSeqST;
      } /*If: finding direction with kmers*/

      seqToIndex_alnSet(qryStackST.seqStr);

      /**************************************************\
      * Main Sec03 Sub02:
      *   - align and print out alignment
      \**************************************************/

      scoreSL =
         memwater(
            &qryStackST,
            refSTPtr,
            &refStartSL,
            &refEndSL,
            &qryStartSL,
            &qryEndSL,
            &setStackST
         );

      fprintf(
         outFILE,
         "%s\t%s\t%0.2f",
         qryStackST.idStr,
         refSTPtr->idStr,
         (float) scoreSL / (float) def_scoreAdj_alnDefs
      ); /*non-coodinates*/

      if(! kmerScanBl)
         goto noScan_main_sec03_sub02;

      if(numKmersSI < 0)
      { /*If: read is reverse complement*/
         fprintf(
            outFILE,
            "\tR\t%lu\t%lu\t%lu\t%lu%s",
            refSTPtr->seqLenSL - refEndSL,
            refSTPtr->seqLenSL - refStartSL,
            qryStackST.seqLenSL - qryEndSL,
            qryStackST.seqLenSL - qryStartSL,
            str_endLine
         ); /*need to convert coordinates*/
         /*seqLenSL is index 1, while xEndUL and xStartUL
         `   are index 0. So, result is index 1
         */
      } /*If: read is reverse complement*/

      else
      { /*Else: read is forward*/
         noScan_main_sec03_sub02:;

         fprintf(
            outFILE,
            "\tF\t%lu\t%lu\t%lu\t%lu%s",
            refStartSL + 1,
            refEndSL + 1,
            qryStartSL + 1,
            qryEndSL + 1,
            str_endLine
         );
      } /*Else: read is forward*/

      /**************************************************\
      * Main Sec03 Sub03:
      *   - get next query sequence
      \**************************************************/

      if(qryTypeSC == def_fqFile_mainMemwater)
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

      goto memErr_main_sec04_sub02;
   } /*If: memory error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - clean up
   ^   o main sec04 sub01:
   ^     - no error cases
   ^   o main sec04 sub02:
   ^     - memory error cases
   ^   o main sec04 sub03:
   ^     - file error cases
   ^   o main sec04 sub04:
   ^     - clean up (error or no error)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - no error cases
   \*****************************************************/

   errSC = 0;
   goto cleanUp_main_sec04_sub04;

   /*****************************************************\
   * Main Sec04 Sub02:
   *   - memory error cases
   \*****************************************************/

   memErr_main_sec04_sub02:;
   errSC = 1;
   goto cleanUp_main_sec04_sub04;
   
   /*****************************************************\
   * Main Sec04 Sub03:
   *   - file error cases
   \*****************************************************/

   fileErr_main_sec04_sub03:;
   errSC = 2;
   goto cleanUp_main_sec04_sub04;

   /*****************************************************\
   * Main Sec04 Sub04:
   *   - clean up (error or no error)
   \*****************************************************/

   cleanUp_main_sec04_sub04:;

   if(
         seqFILE
      && seqFILE != stdin
      && seqFILE != stdout
   ) fclose(seqFILE);

   seqFILE = 0;

   if(
         outFILE
      && outFILE != stdin
      && outFILE != stdout
   ) fclose(outFILE);

   outFILE = 0;
     
   freeStack_seqST(&qryStackST);
   freeStack_alnSet(&setStackST);
   freeStack_kmerCnt(&kmerStackST);

   free(kmerHeapArySI);
   kmerHeapArySI = 0;

   free(cntHeapArySI);
   cntHeapArySI = 0;

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
: Copyright (c) 2024 jeremyButtler
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
