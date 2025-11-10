/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' SOF: Start Of File
'   - translate sequence(s) in to amino acid sequences
'   o header:
'     - included libraries and defaults
'   o fun01: pversion_tranSeq
'     - prints tranSeq's version numbers
'   o fun02: phelp_tranSeq
'     - prints out tranSeq's help message
'   o fun03: input_tranSeq
'     - gets user input
'   o main:
'     - translate sequence(s) into amino acids
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*--------------------------------------------------------\
| Header:
|  - Includes
\--------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../genLib/base10str.h"
#include "../genLib/ulCp.h"
#include "../genBio/codonFun.h"
#include "../genBio/seqST.h"

/*only .h files*/
#include "../bioTools.h" /*version number*/
#include "../genLib/endLine.h" /*version number*/

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - .c  #include "../genLib/fileFun.h"
!   - .h  #include "../genBio/ntTo2Bit.h"
!   - .h  #include "../genBio/revNtTo2Bit.h"
!   - .h  #include "../genBio/codonTbl.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_ORF1_tranSeq 1
#define def_ORF2_tranSeq 0
#define def_ORF3_tranSeq 0
#define def_ORF4_tranSeq 0
#define def_ORF5_tranSeq 0
#define def_ORF6_tranSeq 0

/*--------------------------------------------------------\
| Fun01: pversion_tranSeq
|  - prints tranSeq's version numbers
| Input:
|  - outFILE:
|    o pointer to FILE to print to
| Output:
|  - Prints:
|    o version number to outFILE
\--------------------------------------------------------*/
void
pversion_tranSeq(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "tranSeq from bioTools version: %i-%02i-%02i%s",
      def_year_bioTools,
      def_month_bioTools,
      def_day_bioTools,
      str_endLine
   );
} /*pversion_tranSeq*/

/*--------------------------------------------------------\
| Fun02: phelp_tranSeq
|  - prints out tranSeq's help message
| Input:
|  - outFILE:
|    o pointer to FILE to print to
| Output:
|  - Prints:
|    o help message to outFILE
\--------------------------------------------------------*/
void
phelp_tranSeq(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - fun01 sec01:
   '     - usage line
   '   o fun01 sec01:
   '     - input
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - usage line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "tranSeq -fa sequence.fasta <options ...>%s",
      str_endLine
   );

   fprintf(
     (FILE *) outFILE,
     "  - translates all sequences in a fasta file%s",
     str_endLine
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - input
   ^   o fun01 sec02 sub01:
   ^     - input block
   ^   o fun01 sec02 sub02:
   ^     - file IO
   ^   o fun01 sec02 sub03:
   ^     - coordinate input
   ^   o fun01 sec02 sub04:
   ^     - ofr (open reading frame) selection
   ^   o fun01 sec02 sub05:
   ^     - help message and version number
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *   - input block
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "Input:%s",
      str_endLine
   );

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *   - file IO
   *   o fun01 sec02 sub02 cat01:
   *     - input fasta file
	*   o fun01 sec02 sub02 cat02:
   *     - output file
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub02 Cat01:
   +   - input fasta file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      (FILE *) outFILE,
      "  -fa sequence.fasta: [Required]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o fasta file with sequence(s) to translate%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o use \"-fa -\" for stdin input%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub02 Cat02:
   +   - output file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      (FILE *) outFILE,
      "  -out output.fasta: [Option; stdout]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o fasta file to save translated sequences to%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o use \"-out -\" for stdout input%s",
      str_endLine
   );

   /*****************************************************\
   * Fun01 Sec02 Sub03:
   *   - coordinate input
   *   o fun01 sec02 sub03 cat01:
   *     - starting coordinate input
   *   o fun01 sec02 sub03 cat02:
   *     - ending coordinate input
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub03 Cat01:
   +   - starting coordinate input
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      (FILE *) outFILE,
      "  -start 0: [Optional 0]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o first base to translate in sequence(s)%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub03 Cat02:
   +   - ending coordinate input
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   fprintf(
      (FILE *) outFILE,
      "  -end 0: [Optional 0]%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o last base to translate in sequence(s)%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o use \"-end 0\" for all bases%s",
      str_endLine
   );

   /*****************************************************\
   * Fun01 Sec02 Sub04:
   *   - ofr (open reading frame) selection
   *   o fun01 sec02 sub03 cat01:
   *     - orf 1
   *   o fun01 sec02 sub03 cat02:
   *     - orf 2
   *   o fun01 sec02 sub03 cat03:
   *     - orf 3
   *   o fun01 sec02 sub03 cat04:
   *     - orf 4
   *   o fun01 sec02 sub03 cat05:
   *     - orf 5
   *   o fun01 sec02 sub03 cat06:
   *     - orf 6
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub04 Cat01:
   +   - orf 1
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_ORF1_tranSeq)
      fprintf(
         (FILE *) outFILE,
         "  -orf-1: [Optional Yes]%s",
         str_endLine
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -orf-1: [Optional No]%s",
      str_endLine
      );

   fprintf(
      (FILE *) outFILE,
      "   o translate first reading frame%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o disable with \"-no-orf-1\"%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub04 Cat02:
   +   - orf 2
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_ORF2_tranSeq)
      fprintf(
         (FILE *) outFILE,
         "  -orf-2: [Optional Yes]%s",
      str_endLine
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -orf-2: [Optional No]%s",
      str_endLine
      );

   fprintf(
      (FILE *) outFILE,
      "   o translate second reading frame%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o disable with \"-no-orf-2\"%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub04 Cat03:
   +   - orf 3
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_ORF3_tranSeq)
      fprintf(
         (FILE *) outFILE,
         "  -orf-3: [Optional Yes]%s",
      str_endLine
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -orf-3: [Optional No]%s",
      str_endLine
      );

   fprintf(
      (FILE *) outFILE,
      "   o translate third reading frame%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o disable with \"-no-orf-3\"%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub04 Cat04:
   +   - orf 4
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_ORF4_tranSeq)
      fprintf(
         (FILE *) outFILE,
         "  -orf-4: [Optional Yes]%s",
      str_endLine
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -orf-4: [Optional No]%s",
      str_endLine
      );

   fprintf(
      (FILE *) outFILE,
      "   o translate fourth reading frame (reverse)%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o disable with \"-no-orf-4\"%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub04 Cat05:
   +   - orf 5
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_ORF5_tranSeq)
      fprintf(
         (FILE *) outFILE,
         "  -orf-5: [Optional Yes]%s",
      str_endLine
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -orf-5: [Optional No]%s",
      str_endLine
      );

   fprintf(
      (FILE *) outFILE,
      "   o translate fith reading frame (reverse)%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o disable with \"-no-orf-5\"%s",
      str_endLine
   );

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub04 Cat06:
   +   - orf 6
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(def_ORF6_tranSeq)
      fprintf(
         (FILE *) outFILE,
         "  -orf-6: [Optional Yes]%s",
      str_endLine
      );
   else
      fprintf(
         (FILE *) outFILE,
         "  -orf-6: [Optional No]%s",
      str_endLine
      );

   fprintf(
      (FILE *) outFILE,
      "   o translate six reading frame (reverse)%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "   o disable with \"-no-orf-6\"%s",
      str_endLine
   );

   /*****************************************************\
   * Fun01 Sec02 Sub05:
   *   - help message and version number
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
   ^ Fun01 Sec03:
   ^   - output
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   fprintf(
      (FILE *) outFILE,
      "Output:%s",
      str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  - prints amino acid sequence(s) to stdout%s",
      str_endLine
   );
} /*phelp_tranSeq*/

/*-------------------------------------------------------\
| Fun03: input_tranSeq
|   - gets user input
| Input:
|   - numArgsSI:
|     o number arguments user input
|   - argAryStr:
|     o c-string array with user input
|   - faStrPtr:
|     o c-string pointer to update with fasta file path
|   - outStrPtr:
|     o c-string pointer to update with output file path
|   - startSLPtr:
|     o signed long pointer to hold starting coordinate
|   - endSLPtr:
|     o signed long pointer to hold endig coordinate
|   - orf1BlPtr:
|     o signed char pointer
|       * set to 1 if printing ORF1 results
|       * set to 0 if note printing ORF1 results
|   - orf2BlPtr:
|     o signed char pointer
|       * set to 1 if printing ORF2 results
|       * set to 0 if note printing ORF2 results
|   - orf3BlPtr:
|     o signed char pointer
|       * set to 1 if printing ORF3 results
|       * set to 0 if note printing ORF3 results
|   - orf4BlPtr:
|     o signed char pointer
|       * set to 1 if printing ORF4 results
|       * set to 0 if note printing ORF4 results
|   - orf5BlPtr:
|     o signed char pointer
|       * set to 1 if printing ORF5 results
|       * set to 0 if note printing ORF5 results
|   - orf6BlPtr:
|     o signed char pointer
|       * set to 1 if printing ORF6 results
|       * set to 0 if note printing ORF6 results
| Output:
|   - Modifies:
|     o all input variables to have user input
|   - Prints:
|     o if had error; error message to stderr
|   - Returns:
|     o 0 for no errors
|     o 1 if printed help message
|     o 2 if had error
\-------------------------------------------------------*/
signed char
input_tranSeq(
   int numArgsSI,           /*number arguments input*/
   char *argAryStr[],       /*arguments user input*/
   signed char **faStrPtr,  /*gets fasta file*/
   signed char **outStrPtr, /*gets output file*/
   signed long *startSLPtr, /*gets starting position*/
   signed long *endSLPtr,   /*gets ending position*/
   signed char *orf1BlPtr,  /*boolean for ORF 1*/ 
   signed char *orf2BlPtr,  /*boolean for ORF 2*/ 
   signed char *orf3BlPtr,  /*boolean for ORF 3*/ 
   signed char *orf4BlPtr,  /*boolean for ORF 4*/ 
   signed char *orf5BlPtr,  /*boolean for ORF 5*/ 
   signed char *orf6BlPtr   /*boolean for ORF 6*/ 
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC:
   '   - gets user input
   '   o fun03 sec01:
   '     - variable declartions + input check
   '   o fun03 sec02:
   '     - get input
   '   o fun03 sec03:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec01:
   ^   - variable declartions + input check
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int siArg = 1;
   signed char errSC = 0;
   signed char *tmpStr = 0;

   if(numArgsSI < 1)
      goto phelp_fun03_sec03;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec02:
   ^   - get input
   ^   o fun03 sec02 sub01:
   ^     - file input/output + start loop
   ^   o fun03 sec02 sub02:
   ^     - get starting and ending coordinate
   ^   o fun03 sec02 sub03:
   ^     - get open reading frame(s) (orf) to translate
   ^   o fun03 sec02 sub04:
   ^     - help message requests
   ^   o fun03 sec02 sub05:
   ^     - version number requests
   ^   o fun03 sec02 sub06:
   ^     - invalid input
   ^   o fun03 sec02 sub07:
   ^     - move to next argument
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun03 Sec02 Sub01:
   *   - file input/output + start loop
   *   o fun03 sec02 sub01 cat01:
   *     - input fasta file + start loop
   *   o fun03 sec02 sub01 cat02:
   *     - output file
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun03 Sec02 Sub01 Cat01:
   +   - input fasta file + start loop
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   while(siArg < numArgsSI)
   { /*Loop: get user input*/

      if(
         ! eqlNull_ulCp(
            (signed char *) "-fa",
            (signed char *) argAryStr[siArg]
         )
      ){ /*If: fasta file input*/
         ++siArg;
         *faStrPtr = (signed char *) argAryStr[siArg];
      }  /*If: fasta file input*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub01 Cat02:
      +   - output file
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-out",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: output file input*/
         ++siArg;
         *outStrPtr = (signed char *) argAryStr[siArg];
      }  /*Else If: output file input*/

      /**************************************************\
      * Fun03 Sec02 Sub02:
      *   - get starting and ending coordinate
      *   o fun03 sec02 sub02 cat01:
      *     - get starting coordinate
      *   o fun03 sec02 sub02 cat02:
      *     - get ending coordinate
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub02 Cat01:
      +   - get starting coordinate
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-start",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: start coordinate input*/
         ++siArg;
         tmpStr = (signed char *) argAryStr[siArg];

         tmpStr +=
            strToSL_base10str(
               tmpStr,
               startSLPtr
            ); /*convert string to number*/

         if(*tmpStr != '\0')
         { /*If: error in conversion*/
            fprintf(
               stderr,
               "-start %s is non-numeric or to large%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec03;
         } /*If: error in conversion*/
      } /*Else If: start coordinate input*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub02 Cat02:
      +   - get ending coordinate
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-end",
            (signed char *) argAryStr[siArg]
         )
      ){ /*Else If: end coordinate input*/
         ++siArg;
         tmpStr = (signed char *) argAryStr[siArg];

         tmpStr +=
            strToSL_base10str(
               tmpStr,
               endSLPtr
            ); /*convert string to number*/

         if(*tmpStr != '\0')
         { /*If: error in conversion*/
            fprintf(
               stderr,
               "-end %s is non-numeric or to large%s",
               argAryStr[siArg],
               str_endLine
            );

            goto err_fun03_sec03;
         } /*If: error in conversion*/
      } /*Else If: end coordinate input*/

      /**************************************************\
      * Fun03 Sec02 Sub03:
      *   - get open reading frame(s) (orf) to translate
      *   o fun03 sec02 sub03 cat01:
      *     - orf 1
      *   o fun03 sec02 sub03 cat02:
      *     - orf 2
      *   o fun03 sec02 sub03 cat03:
      *     - orf 3
      *   o fun03 sec02 sub03 cat04:
      *     - orf 4
      *   o fun03 sec02 sub03 cat05:
      *     - orf 5
      *   o fun03 sec02 sub03 cat06:
      *     - orf 6
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub03 Cat01:
      +   - orf 1
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-orf-1",
            (signed char *) argAryStr[siArg]
         )
      ) *orf1BlPtr = 1;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-orf-1",
            (signed char *) argAryStr[siArg]
         )
      ) *orf1BlPtr = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub03 Cat02:
      +   - orf 2
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-orf-2",
            (signed char *) argAryStr[siArg]
         )
      ) *orf2BlPtr = 1;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-orf-2",
            (signed char *) argAryStr[siArg]
         )
      ) *orf2BlPtr = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub03 Cat03:
      +   - orf 3
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-orf-3",
            (signed char *) argAryStr[siArg]
         )
      ) *orf3BlPtr = 3;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-orf-3",
            (signed char *) argAryStr[siArg]
         )
      ) *orf3BlPtr = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub03 Cat04:
      +   - orf 4
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-orf-4",
            (signed char *) argAryStr[siArg]
         )
      ) *orf4BlPtr = 1;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-orf-4",
            (signed char *) argAryStr[siArg]
         )
      ) *orf4BlPtr = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub03 Cat05:
      +   - orf 5
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-orf-5",
            (signed char *) argAryStr[siArg]
         )
      ) *orf5BlPtr = 1;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-orf-5",
            (signed char *) argAryStr[siArg]
         )
      ) *orf5BlPtr = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun03 Sec02 Sub03 Cat06:
      +   - orf 6
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-orf-6",
            (signed char *) argAryStr[siArg]
         )
      ) *orf6BlPtr = 1;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-no-orf-6",
            (signed char *) argAryStr[siArg]
         )
      ) *orf6BlPtr = 0;

      /**************************************************\
      * Fun03 Sec02 Sub04:
      *   - help message requests
      \**************************************************/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-h",
            (signed char *) argAryStr[siArg]
         )
      ) goto phelp_fun03_sec03;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "--h",
            (signed char *) argAryStr[siArg]
         )
      ) goto phelp_fun03_sec03;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "help",
            (signed char *) argAryStr[siArg]
         )
      ) goto phelp_fun03_sec03;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-help",
            (signed char *) argAryStr[siArg]
         )
      ) goto phelp_fun03_sec03;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "--help",
            (signed char *) argAryStr[siArg]
         )
      ) goto phelp_fun03_sec03;


      /**************************************************\
      * Fun03 Sec02 Sub05:
      *   - version number requests
      \**************************************************/

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-v",
            (signed char *) argAryStr[siArg]
         )
      ) goto pversion_fun03_sec03;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "--v",
            (signed char *) argAryStr[siArg]
         )
      ) goto pversion_fun03_sec03;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "version",
            (signed char *) argAryStr[siArg]
         )
      ) goto pversion_fun03_sec03;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "-version",
            (signed char *) argAryStr[siArg]
         )
      ) goto pversion_fun03_sec03;

      else if(
         ! eqlNull_ulCp(
            (signed char *) "--version",
            (signed char *) argAryStr[siArg]
         )
      ) goto pversion_fun03_sec03;

      /**************************************************\
      * Fun03 Sec02 Sub06:
      *   - invalid input
      \**************************************************/

      else
      { /*Else: unrecognized input*/
         fprintf(
            stderr,
            "%s is not recongnized%s",
            argAryStr[siArg],
            str_endLine
         );

         goto err_fun03_sec03;
      } /*Else: unrecognized input*/

      /**************************************************\
      * Fun03 Sec02 Sub07:
      *   - move to next argument
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun03 Sec03:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun03_sec03;

   phelp_fun03_sec03:;
      errSC = 1;
      phelp_tranSeq(stdout);
      goto ret_fun03_sec03;
   pversion_fun03_sec03:;
      errSC = 1;
      pversion_tranSeq(stdout);
      goto ret_fun03_sec03;
   err_fun03_sec03:;
      errSC = 2;
      goto ret_fun03_sec03;
   ret_fun03_sec03:;
      return errSC;
} /*input_tranSeq*/

/*--------------------------------------------------------\
| Main:
|   - translates sequence(s) into amino acids
| Input:
|   - numArgsI:
|     o Number of arguments the user input
|   - argsStr:
|     o Arguments and parameters the user input
| Output:
|   - Prints:
|     o The translated sequences to stdout
\--------------------------------------------------------*/
int main(
   int numArgSI,    /*number of arguments the user input*/
   char *argAryStr[]/*Arguments/parameters input by user*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main: main
   '   - translates sequence(s) into amino acids
   '   o main sec01:
   '     - Variable declerations
   '   o main sec02:
   '     - Check user input and initialize variables
   '   o main sec03:
   '     - Get and reverse complement fasta entries
   '   o main sec04:
   '     - Clean up and exit
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^   - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *faFileStr = 0;
   signed char *outFileStr = 0;
   signed long startSL = 0;
   signed long endSL = 0;

   signed long tmpStartSL = 0;
   signed long tmpEndSL = 0;
   signed long revEndSL = 0;

   signed char orf1Bl = def_ORF1_tranSeq;
   signed char orf2Bl = def_ORF2_tranSeq;
   signed char orf3Bl = def_ORF3_tranSeq;
   signed char orf4Bl = def_ORF4_tranSeq;
   signed char orf5Bl = def_ORF5_tranSeq;
   signed char orf6Bl = def_ORF6_tranSeq;

   
   signed char errSC = 0;
   signed long numSeqSL = 0;

   signed char *aaHeapStr = 0;
   signed long lenAaSL = 0;
   signed long sizeAaSL = 0;

   FILE *seqFILE = 0;
   FILE *outFILE = 0;

   struct seqST seqStackST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^   - Check user input and initialize variables
   ^   o main sec02 sub01:
   ^     - initialize structs
   ^   o main sec02 sub02:
   ^     - get user input
   ^   o main sec02 sub02:
   ^     - Initialize variables
   ^   o main sec02 sub03:
   ^     - open sequence file
	^   o main sec02 sub04:
   ^     - open output file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Main Sec02 Sub01:
   *   - initialize structs
   \******************************************************/

   init_seqST(&seqStackST);

   /******************************************************\
   * Main Sec02 Sub02:
   *   - get user input
   \******************************************************/

   errSC =
      input_tranSeq(
         numArgSI,
         argAryStr,
         &faFileStr,
         &outFileStr,
         &startSL,
         &endSL,
         &orf1Bl,
         &orf2Bl,
         &orf3Bl,
         &orf4Bl,
         &orf5Bl,
         &orf6Bl
   );
    

   if(errSC)
   { /*If: had error*/
      --errSC; /*convert help/version to no error (0)*/
      goto ret_main_sec0x;
   } /*If: had error*/

   /******************************************************\
   * Main Sec02 Sub03:
   *   - open sequence file
   \******************************************************/

   if(
         ! faFileStr
      || faFileStr[0] == '-'
   ) seqFILE = stdin;

   else
   { /*Else: user input fasta file*/
      seqFILE =
         fopen(
            (char *) faFileStr,
            "r"
         );

      if(! seqFILE)
      { /*If: could no open fasta file*/
         fprintf(
            stderr,
            "could not open -fa %s%s",
            faFileStr,
            str_endLine
         );

         goto err_main_sec0x;
      } /*If: could not open fasta file*/
   } /*Else: user input fasta file*/

   if(! faFileStr)
      faFileStr = (signed char *) "-";

   /******************************************************\
   * Main Sec02 Sub04:
   *   - open output file
   \******************************************************/

   if(
         ! outFileStr
      || outFileStr[0] == '-'
   ) outFILE = stdout;

   else
   { /*Else: user input output file*/
      outFILE =
         fopen(
            (char *) outFileStr,
            "w"
         );

      if(! seqFILE)
      { /*If: could no open output file*/
         fprintf(
            stderr,
            "could not open -out %s%s",
            outFileStr,
            str_endLine
         );

         goto err_main_sec0x;
      } /*If: could not open output file*/
   } /*Else: user input output file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^  - translate all sequences in the fasta file
   ^  o main sec03 sub01:
   ^    - read in first sequence + start loop
   ^  o main sec03 sub02:
   ^    - start loop (memory allocation part of loop)
   ^  o main sec03 sub03:
   ^    - translate reading frame 1
   ^  o main sec03 sub04:
   ^    - translate reading frame 2
   ^  o main sec03 sub05:
   ^    - translate reading frame 3
   ^  o main sec03 sub06:
   ^    - reverse complement the sequence for ORFs 4 to 6
   ^  o main sec03 sub07:
   ^    - translate reading frame 4
   ^  o main sec03 sub08:
   ^    - translate reading frame 5
   ^  o main sec03 sub09:
   ^    - translate reading frame 6
   ^  o main sec03 sub10:
   ^    - read in the next sequence to translate
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /******************************************************\
   * Main Sec03 Sub01:
   *  - read in first sequence
   \******************************************************/

   errSC = getFa_seqST(seqFILE, &seqStackST);

   aaHeapStr =
      malloc(
         (seqStackST.seqLenSL + 9) * sizeof(signed char)
      );
   if(! aaHeapStr)
      goto memErr_main_sec0x;
   sizeAaSL = seqStackST.seqLenSL;

   if(errSC == def_EOF_seqST)
   { /*If: read end of file*/
      if(seqStackST.seqStr[0] != '\0')
         errSC = 0;
      else
      { /*Else: nothing in file*/
         fprintf(
            stderr,
            "-fa %s had no sequences%s",
            faFileStr,
            str_endLine
         );
      } /*Else: nothing in file*/

      goto err_main_sec0x;
   } /*If: read end of file*/

   else if(errSC)
      goto readErr_main_sec0x;


   /******************************************************\
   * Main Sec03 Sub02:
   *  - start loop (memory allocation part of loop)
   \******************************************************/

   while(! errSC)
   { /*Loop: translate all sequences in fasta file*/
      ++numSeqSL;

      if(sizeAaSL < (signed long) seqStackST.seqLenSL)
      { /*If: need more memory*/
         free(aaHeapStr);
         aaHeapStr = 0;

         aaHeapStr =
            malloc(
                 (seqStackST.seqLenSL + 9)
               * sizeof(signed char)
            );
         if(! aaHeapStr)
            goto memErr_main_sec0x;
         sizeAaSL = seqStackST.seqLenSL;
      } /*If: need more memory*/

      if(seqStackST.seqLenSL < startSL)
         goto nextSeq_main_sec03_sub09;
         /*can not translate*/

      tmpStartSL = startSL;
      if(seqStackST.seqLenSL < endSL)
         tmpEndSL = seqStackST.seqLenSL;
      else if(! endSL)
         tmpEndSL = seqStackST.seqLenSL;
      else
         tmpEndSL = endSL;

      /***************************************************\
      * Main Sec03 Sub03:
      *  - translate reading frame 1
      \***************************************************/

      if(orf1Bl)
      { /*If: translating orf 1*/
         lenAaSL =
            seqToAA_codonFun(
               seqStackST.seqStr,
               aaHeapStr,
               tmpStartSL,
               tmpEndSL
            ); /*Get the aa sequence for ORF 1*/
      
         if(lenAaSL == def_unkownNt_codonFun)
            goto invalidNt_main_sec0x;
         /*else no error or spare bases at end*/

         fprintf(
            (FILE *) outFILE,
            ">%s-ORF1\tstart=%li\tend=%li%s%s%s",
            seqStackST.idStr,
            tmpStartSL,
            tmpEndSL,
            str_endLine,
            aaHeapStr,
            str_endLine
         );
      } /*If: translating orf 1*/

      /***************************************************\
      * Main Sec03 Sub04:
      *  - Translate reading frame 2
      \***************************************************/

      if(orf2Bl)
      { /*If: I am working on reading frame 2*/
         lenAaSL =
            seqToAA_codonFun(
               seqStackST.seqStr,
               aaHeapStr,
               tmpStartSL + 1,
               tmpEndSL
            );
      
         if(lenAaSL == def_unkownNt_codonFun)
            goto invalidNt_main_sec0x;
         /*else no error or spare bases at end*/

         fprintf(
            (FILE *) outFILE,
            ">%s-ORF2\tstart=%li\tend=%li%s%s%s",
            seqStackST.idStr,
            tmpStartSL + 1,
            tmpEndSL,
            str_endLine,
            aaHeapStr,
            str_endLine
         );
      } /*If: I am working on reading frame 2*/

      /***************************************************\
      * Main Sec03 Sub05:
      *  - Translate reading frame 3
      \***************************************************/

      if(orf3Bl)
      { /*If: I am working on reading frame 3*/
         lenAaSL =
            seqToAA_codonFun(
               seqStackST.seqStr,
               aaHeapStr,
               tmpStartSL + 2,
               tmpEndSL
            );
      
         if(lenAaSL == def_unkownNt_codonFun)
            goto invalidNt_main_sec0x;
         /*else no error or spare bases at end*/

         fprintf(
            (FILE *) outFILE,
            ">%s-ORF3\tstart=%li\tend=%li%s%s%s",
            seqStackST.idStr,
            tmpStartSL + 2,
            tmpEndSL,
            str_endLine,
            aaHeapStr,
            str_endLine
         );
      } /*If: I am working on reading frame 3*/

      /***************************************************\
      * Main Sec03 Sub06:
      *  - Reverse complement the sequence for ORFs 4 to 6
      \***************************************************/

      /*Check if I am printing the reverse reading frames*/
      if(orf4Bl || orf5Bl || orf6Bl)
      { /*If: need to reverse complment the sequence*/
         revComp_seqST(&seqStackST);

         if(! startSL && ! endSL)
            ; /*entire sequence no conversion needed*/

         else if(startSL && endSL)
         { /*Else If: adjust and swap coordinates*/
            /*need to refind stat and end coordinates*/
            tmpStartSL = seqStackST.seqLenSL - tmpEndSL;
            revEndSL = seqStackST.seqLenSL - startSL;
         } /*Else If: adjust and swap coordinates*/

         else if(startSL)
         { /*Else If: only using starting coordinate*/
            revEndSL = seqStackST.seqLenSL - startSL;
            tmpStartSL = 0;
         } /*Else If: only using starting coordinate*/

         else
         { /*Else: only using ending coordinate*/
            tmpStartSL = seqStackST.seqLenSL - tmpEndSL;
            revEndSL = 0;
         } /*Else: only using ending coordinate*/
      }  /*If: need to reverse complment the sequence*/

      else
         goto nextSeq_main_sec03_sub09;

      /***************************************************\
      * Main Sec03 Sub07:
      *  - Translate reading frame 4
      \***************************************************/

      if(orf4Bl)
      { /*If: I am working on reading frame 4*/
         lenAaSL =
            seqToAA_codonFun(
               seqStackST.seqStr,
               aaHeapStr,
               tmpStartSL,
               revEndSL
            );
      
         if(lenAaSL == def_unkownNt_codonFun)
            goto invalidNt_main_sec0x;
         /*else no error or spare bases at end*/

         fprintf(
            (FILE *) outFILE,
            ">%s-ORF4\trev_start=%li\trev_end=%li",
            seqStackST.idStr,
            tmpStartSL,
            revEndSL
         );
         fprintf(
            (FILE *) outFILE,
            "\tfor_start=%li\tfor_end=%li%s%s%s",
            startSL,
            tmpEndSL,
            str_endLine,
            aaHeapStr,
            str_endLine
         );
      } /*If: I am working on reading frame 4*/

      /***************************************************\
      * Main Sec03 Sub08:
      *  - Translate reading frame 5
      \***************************************************/

      if(orf5Bl)
      { /*If: I am working on reading frame 5*/
         lenAaSL =
            seqToAA_codonFun(
               seqStackST.seqStr,
               aaHeapStr,
               tmpStartSL + 1,
               revEndSL
            );
      
         if(lenAaSL == def_unkownNt_codonFun)
            goto invalidNt_main_sec0x;
         /*else no error or spare bases at end*/

         fprintf(
            (FILE *) outFILE,
            ">%s-ORF5\trev_start=%li\trev_end=%li",
            seqStackST.idStr,
            tmpStartSL + 1,
            revEndSL
         );

         fprintf(
            (FILE *) outFILE,
            "\tfor_start=%li\tfor_end=%li%s%s%s",
            startSL,
            tmpEndSL - 1,
            str_endLine,
            aaHeapStr,
            str_endLine
         );
      } /*If: I am working on reading frame 5*/

      /***************************************************\
      * Main Sec03 Sub09:
      *  - Translate reading frame 6
      \***************************************************/

      if(orf6Bl)
      { /*If: I am working on reading frame 6*/
         lenAaSL =
            seqToAA_codonFun(
               seqStackST.seqStr,
               aaHeapStr,
               tmpStartSL + 2,
               revEndSL
            );
      
         if(lenAaSL == def_unkownNt_codonFun)
            goto invalidNt_main_sec0x;
         /*else no error or spare bases at end*/

         fprintf(
            (FILE *) outFILE,
            ">%s-ORF6\trev_start=%li\trev_end=%li",
            seqStackST.idStr,
            tmpStartSL + 2,
            revEndSL
         );
         fprintf(
            (FILE *) outFILE,
            "\tfor_start=%li\tfor_end=%li%s%s%s",
            startSL,
            tmpEndSL - 2,
            str_endLine,
            aaHeapStr,
            str_endLine
         );
      } /*If: I am working on reading frame 5*/

      /***************************************************\
      * Main Sec03 Sub10:
      *  - Read in the next sequence to translate
      \***************************************************/

      nextSeq_main_sec03_sub09:;
         blank_seqST(&seqStackST);
         errSC = getFa_seqST(seqFILE, &seqStackST);
   } /*Loop: translate all sequences in fasta file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^  - Clean up and exit
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(errSC != def_EOF_seqST)
      goto readErr_main_sec0x;

   errSC = 0;
   goto ret_main_sec0x;

   memErr_main_sec0x:;
      fprintf(
         stdout,
         "memory error on sequence %li in -fa %s%s",
         numSeqSL,
         faFileStr,
         str_endLine
      ); /*Print out the error*/

      goto err_main_sec0x;

   readErr_main_sec0x:;
      fprintf(
         stdout,
         "sequence %li in -fa %s is an invalid entry%s",
         numSeqSL,
         faFileStr,
         str_endLine
      ); /*Print out the error*/

      goto err_main_sec0x;

   invalidNt_main_sec0x:;
      fprintf(
         stdout,
         "sequence %li in -fa %s has an unknown base%s",
         numSeqSL,
         faFileStr,
         str_endLine
      ); /*Print out the error*/

      goto err_main_sec0x;

   err_main_sec0x:;
      errSC = 1;
      goto ret_main_sec0x;

   ret_main_sec0x:;
      if(aaHeapStr)
         free(aaHeapStr);
      aaHeapStr = 0;

      freeStack_seqST(&seqStackST);

      if(
            seqFILE
         && seqFILE != stdin
         && seqFILE != stdout
         && seqFILE != stderr
      ) fclose(seqFILE);
      seqFILE = 0;

      if(
            outFILE
         && outFILE != stdin
         && outFILE != stdout
         && outFILE != stderr
      ) fclose(outFILE);
      outFILE = 0;

      return errSC;
} /*main*/
