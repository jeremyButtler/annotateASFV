/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' mainSamToAln SOF: Start Of File
'   - driver function to convert sam file to aln
'   o header:
'     - included libraries and defined variables
'   o fun01: pversion_mainSamToAln
'     - prints versoin number for samToAln
'   o fun02: phelp_mainSamToAln
'     - prints help message for samToAln
'   o fun03: checkIfHelp_mainSamToAln
'     - checks if input was the help message
'   o fun04: checkIfVersion_mainSamToAln
'     - checks if input was the version number request
'   o fun05: input_mainSamToAln
'     - gets user input
'   o main:
'     - driver function for samToAln
'   o license:
'     - Licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and defined variables
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#include "../genLib/base10str.h"
#include "../genLib/ulCp.h"
#include "../genLib/charCp.h"

#include "../genBio/seqST.h"
#include "../genBio/samEntry.h"

#include "../genAln/alnSet.h"
#include "../genAln/samToAln.h"

/*no .c files*/
#include "../genLib/endLine.h"
#include "../bioTools.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden dependices
!   - .c  #include "../genLib/numToStr.h"
!   - .c  #include "../genLib/fileFun.h"
!   - .h  #include "../genBio/ntTo5Bit.h"
!   - .h  #include "../genAln/alnDefs.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_pmask_mainSamToAln 0    /*print masked bases*/
#define def_ppos_mainSamToAln 1     /*print positions*/
#define def_lineWrap_mainSamToAln 58 /*line wrap*/
#define def_pFullId_mainSamToAln 1  /*print full ref id*/

/*-------------------------------------------------------\
| Fun01: pversion_mainSamToAln
|   - prints versoin number for samToAln
| Input:
|   - outFILE:
|     o file to print version  number to
| Output:
|   - Prints:
|     o version number to outFILE
\-------------------------------------------------------*/
void
pversion_mainSamToAln(
   void *outFILE
){
   fprintf(
      (FILE *) outFILE,
      "samToAln from bioTools version: %i-%02i-%02i%s",
      def_year_bioTools,
      def_month_bioTools,
      def_day_bioTools,
      str_endLine
   );
} /*pversion_mainSamToAln*/

/*-------------------------------------------------------\
| Fun02: phelp_mainSamToAln
|   - prints help message for samToAln
| Input:
|   - outFILE:
|     o file to print help message to
| Output:
|   - Prints:
|     o help message to outFILE
\-------------------------------------------------------*/
void
phelp_mainSamToAln(
   void *outFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC:
   '   - prints help message for samToAln
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
     "samToAln -ref ref.fa -sam reads.sam [options...]%s",
     str_endLine
   );

   fprintf(
      (FILE *) outFILE,
      "  - converts sam file to expand cigar aligment%s",
      str_endLine
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun02 Sec02:
   ^   - print input block
   ^   o fun02 sec02 sub01:
   ^     - input block
   ^   o fun02 sec02 sub02:
   ^     - file input arguments
   ^   o fun02 sec02 sub03:
   ^     - output arguments
   ^   o fun02 sec02 sub04:
   ^     - print help/version number commands
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun02 Sec02 Sub01:
   *   - input block
   \*****************************************************/

   fprintf(
      (FILE *) outFILE,
      "Input:%s",
      str_endLine
   );

  /*****************************************************\
  * Fun02 Sec02 Sub02:
  *   - file input arguments
  \*****************************************************/

   fprintf(
       (FILE *) outFILE,
       "  -ref ref.fa: [Required]%s",
       str_endLine              
    );

   fprintf(
       (FILE *) outFILE,
       "    o fasta with reference sequence used for%s",
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "      mapping%s",
       str_endLine
    );


   /*sam file*/
   fprintf(
       (FILE *) outFILE,
       "  -sam reads.sam: [Required; stdin]%s",
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o sam file with mapped reads%s",
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o use \"-sam -\" for stdin input%s",
       str_endLine
    );


   /*match matrix*/
   fprintf(
       (FILE *) outFILE,
       "  -match matrix.txt: [Optional]%s",
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o match matrix for alignment%s",
       str_endLine
    );

 /*****************************************************\
  * Fun02 Sec02 Sub03:
  *   - output arguments
  \*****************************************************/

   if(def_pFullId_mainSamToAln)
     fprintf(
         (FILE *) outFILE,
         "  -full-id: [Optional; Yes]%s",
         str_endLine
      );
   else
     fprintf(
         (FILE *) outFILE,
         "  -full-id: [Optional; No]%s",
         str_endLine
      );

   fprintf(
      (FILE *) outFILE,
      "    o print refernce id past first space or tab%s",
      str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o use \"-short-id\" to disable%s",
       str_endLine
    );


   fprintf(
       (FILE *) outFILE,
       "  -wrap %i: [Optional]%s",
       def_lineWrap_mainSamToAln,
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o maximum line length (0 is no limit)%s",
       str_endLine
    );


   /*keep/remove masking*/
   if(def_pmask_mainSamToAln)
     fprintf(
         (FILE *) outFILE,
         "  -mask: [Optional; Yes]%s",
         str_endLine
      );
   else
     fprintf(
         (FILE *) outFILE,
         "  -mask: [Optional; No]%s",
         str_endLine
      );

   fprintf(
       (FILE *) outFILE,
       "    o print unaligned reference and query%s",
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o disable with \"-no-mask\"%s",
       str_endLine
    );

   /*print sequence positions*/
   if(def_ppos_mainSamToAln)
     fprintf(
         (FILE *) outFILE,
         "  -pos: [Optional; Yes]%s",
         str_endLine
      );
   else
     fprintf(
         (FILE *) outFILE,
         "  -pos: [Optional; No]%s",
         str_endLine
      );

   fprintf(
       (FILE *) outFILE,
       "    o print out line start and end positions%s",
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o disable with \"-no-pos\"%s",
       str_endLine
    );


   /*output file*/
   fprintf(
       (FILE *) outFILE,
       "  -out out.file: [Optional; stdout]%s",
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o file to save alignmet to%s",
       str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "    o use \"-out -\" for stdout%s",
       str_endLine
    );


  /*****************************************************\
  * Fun02 Sec02 Sub04:
  *   - print help/version number commands
  \*****************************************************/

  fprintf(
      (FILE *) outFILE,
      "  -h: print this help message to stdout and exit%s",
      str_endLine
    );

   fprintf(
       (FILE *) outFILE,
       "  -v: print version number to stdout and exit%s",
       str_endLine
    );

   fprintf(
      (FILE *) outFILE,
      "  -p-match: print default match matrix to stdout%s",
      str_endLine
    );
   fprintf(
       (FILE *) outFILE,
       "    o only used when have \'M\' in cigar%s",
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
      "  - prints alignment to stdout or -out file.aln%s",
      str_endLine
    );
} /*phelp_mainSamToAln*/

/*-------------------------------------------------------\
| Fun03: checkIfHelp_mainSamToAln
|   - checks if input was the help message
| Input:
|   - argStr:
|     o c-string with argument to check for help message
| Output:
|   - Returns:
|     o 1 if is the help message
|     o 0 if not the help message
\-------------------------------------------------------*/
signed char
checkIfHelp_mainSamToAln(
   signed char *argStr
){
   if(! eql_charCp((signed char *) "-h", argStr, '\0'))
     return 1;

   if(! eql_charCp((signed char *) "--h", argStr, '\0'))
     return 1;

   if(! eql_charCp((signed char *) "help", argStr, '\0'))
     return 1;

   if(! eql_charCp((signed char *) "-help", argStr, '\0'))
     return 1;

   if(! eql_charCp((signed char *) "--help",argStr,'\0'))
     return 1;

   return 0;
} /*checkIfHelp_mainSamToAln*/

/*-------------------------------------------------------\
| Fun04: checkIfVersion_mainSamToAln
|   - checks if input was the version number request
| Input:
|   - argStr:
|     o c-string with argument to check for version number
| Output:
|   - Returns:
|     o 1 if is version number
|     o 0 if not version number
\-------------------------------------------------------*/
signed char
checkIfVersion_mainSamToAln(
   signed char *argStr
){
    if(! eql_charCp((signed char *) "-v", argStr, '\0'))
      return 1;

    if(! eql_charCp((signed char *) "--v", argStr, '\0'))
      return 1;

    if(
       ! eql_charCp((signed char *) "version",argStr,'\0')
    ) return 1;

    if(
       ! eql_charCp(
          (signed char *) "-version",argStr,'\0'
       )
    ) return 1;

    if(
       ! eql_charCp(
          (signed char *) "--version", argStr, '\0'
       )
    ) return 1;

    return 0;
} /*checkIfVersion_mainSamToAln*/

/*-------------------------------------------------------\
| Fun05: input_mainSamToAln
|   - gets user input
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o c-string array with user input
|   - samFileStrPtr:
|     o pointer to c-string to point to the sam file
|   - refFileStrPtr:
|     o pointer to c-string to point to refernece file
|   - outFileStrPtr:
|     o pointer to c-string to point to output file
|   - pFullIdBlPtr:
|     o pointer to sigend char to be set to
|       * 1 if user wanted full reference id
|       * 0 if user wanted short (no white space) id
|   - alnSetSTPtr:
|     o pointer to alnSet structure to hold input
| Output:
|   - Returns:
|     o 0 for no errors
|     o 1 if printed help, version number, or match matrix
|     o 2 for errors
\-------------------------------------------------------*/
signed char
input_mainSamToAln(
      int numArgsSI,
      char *argAryStr[],
      signed char **samFileStrPtr,
      signed char **refFileStrPtr,
      signed char **outFileStrPtr,
      signed char *pFullIdBlPtr,
      struct alnSet *alnSetSTPtr
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC:
   '   - gets user input
   '   o fun05 sec02:
   '     - variable declerations
   '   o fun05 sec02:
   '     - check if user input anything
   '   o fun05 sec03:
   '     - get use input
   '   o fun05 sec03:
   '     - return any errors
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   signed int siArg = 1;
   signed char *tmpStr = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - check if user input anything
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(numArgsSI <= 1)
   { /*If: no arguments input*/
     phelp_mainSamToAln(stdout);
     goto phelp_fun05_sec04;
   } /*If: no arguments input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec03:
   ^   - get use input
   ^   o fun05 sec03 sub01:
   ^     - get file io and start get input loop
   ^   o fun05 sec03 sub02:
   ^     - get output formating
   ^   o fun05 sec03 sub03:
   ^     - see if help/version/match settings print
   ^   o fun05 sec03 sub04:
   ^     - else; unkown input
   ^   o fun05 sec03 sub05:
   ^     - move to next argument
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun05 Sec03 Sub01:
   *   - get file io and start get input loop
   \*****************************************************/

    while(siArg < numArgsSI)
    { /*Loop: get user input*/
      if(
         ! eql_charCp(
           (signed char *) "-ref",
           (signed char *) argAryStr[siArg],
           '\0'
         )
      ){ /*If: reference seqeunce*/
        ++siArg;
        *refFileStrPtr = (signed char *) argAryStr[siArg];
      } /*If: reference seqeunce*/

      else if(
         ! eql_charCp(
            (signed char *) "-sam",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ){ /*Else If: sam file with mapped reads*/
        ++siArg;
        *samFileStrPtr = (signed char *) argAryStr[siArg];
      } /*Else If: sam file with mapped reads*/

      else if(
         ! eql_charCp(
            (signed char *) "-out",
            (signed char *) argAryStr[siArg],
            '\0'
        )
      ){ /*Else If: output file*/
        ++siArg;
        *outFileStrPtr = (signed char *) argAryStr[siArg];
      } /*Else If: output file*/

      else if(
         ! eql_charCp(
            (signed char *) "-full-id",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ) *pFullIdBlPtr = 1;
         /*Else If: user wanted full reference id*/

      else if(
         ! eql_charCp(
            (signed char *) "-short-id",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ) *pFullIdBlPtr = 0;
         /*user wanted short, nowhite space reference id*/

      else if(
         ! eql_charCp(
            (signed char *) "-match",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ){ /*Else If: match matrix file*/
          ++siArg;
          errSC =
              readMatchFile_alnSet(
                 alnSetSTPtr,
                 (signed char *) argAryStr[siArg]
             );

          if(errSC)
          { /*If: had an error*/
             fprintf(
                stderr,
                "error reading -match %s%s",
                argAryStr[siArg],
                str_endLine
             );

             goto err_fun05_sec04;
          } /*If: had an error*/
       } /*Else If: match matrix file*/

      /**************************************************\
      * Fun05 Sec03 Sub02:
      *   - get output formating
      \**************************************************/

      else if(
         ! eql_charCp(
            (signed char *) "-mask",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ) alnSetSTPtr->pFullAlnBl = 1;

      else if(
         ! eql_charCp(
            (signed char *) "-no-mask",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ) alnSetSTPtr->pFullAlnBl = 0;

      else if(
         ! eql_charCp(
            (signed char *) "-pos",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ) alnSetSTPtr->pBasePosBl = 1;

      else if(
         ! eql_charCp(
            (signed char *) "-no-pos",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ) alnSetSTPtr->pBasePosBl = 0;

      else if(
         ! eql_charCp(
            (signed char *) "-wrap",
            (signed char *) argAryStr[siArg],
            '\0'
         )
      ){ /*Else If: line wrap*/
          ++siArg;
          
          tmpStr = (signed char *) argAryStr[siArg];

          tmpStr +=
             strToUI_base10str(
                (signed char *) argAryStr[siArg],
                &alnSetSTPtr->lineWrapUI
             );

          if(*tmpStr != '\0')
          { /*If: non-numeric entry*/
             fprintf(
                stderr,
                "non-numeric or to large -wrap %s%s",
                argAryStr[siArg],
                str_endLine
             );

             goto err_fun05_sec04;
          } /*If: non-numeric entry*/
       } /*Else If: line wrap*/

      /**************************************************\
      * Fun05 Sec03 Sub03:
      *   - see if help/version/match settings print
      \**************************************************/

      else if(
         checkIfHelp_mainSamToAln(
           (signed char *) argAryStr[siArg]
        )
      ){ /*Else If: wanted help message*/
         phelp_mainSamToAln(stdout);
         goto phelp_fun05_sec04;
      } /*Else If: wanted help message*/

      else if(
         checkIfVersion_mainSamToAln(
           (signed char *) argAryStr[siArg]
        )
      ){ /*Else If: wanted version number*/
         pversion_mainSamToAln(stdout);
         goto pversion_fun05_sec04;
      } /*Else If: wanted version numbe*/

      else if(
         ! eql_charCp(
            (signed char *) "-p-match",
            (signed char *) argAryStr[siArg],
            '\0'
           )
       ){ /*Else If: printing default match matrix*/
          pDefMatchMatrix_alnSet(stdout);
          goto pmatchMatrix_fun05_sec04;
       } /*Else If: printing default match matrix*/
      
      /**************************************************\
      * Fun05 Sec03 Sub04:
      *   - else; unkown input
      \**************************************************/

       else
       { /*Else: unkown input*/
          fprintf(
             stderr,
             "%s is not reconginized%s",
             argAryStr[siArg],
             str_endLine
           );

          goto err_fun05_sec04;
       } /*Else: unkown input*/

      /**************************************************\
      * Fun05 Sec03 Sub05:
      *   - move to next argument
      \**************************************************/

      ++siArg;
   } /*Loop: get user input*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec04:
   ^   - return any errors
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto cleanUp_fun05_sec04;

   phelp_fun05_sec04:;
   pmatchMatrix_fun05_sec04:;
   pversion_fun05_sec04:;
   errSC = 1;
   goto cleanUp_fun05_sec04;

   err_fun05_sec04:;
   errSC = 2;
   goto cleanUp_fun05_sec04;
  
   cleanUp_fun05_sec04:;
   return errSC;
} /*input_mainSamToAln*/

/*-------------------------------------------------------\
| Main:
|   - driver function for samToAln
| Input:
|   - numArgsSI:
|     o number of arguments the user input
|   - argAryStr:
|     o c-string array of user input arguments
| Output:
|   - Prints:
|     o aligned sequences to stdout or "-out file.aln"
\-------------------------------------------------------*/
int
main(
   int numArgsSI,
   char *argAryStr[]
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main TOC:
   '   o main sec01:
   '     - variable declarations
   '   o main sec02:
   '     - initialize, get user input, open files
   '   o main sec03:
   '     - print out alignments
   '   o main sec04:
   '     - clean up and return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^    - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char *tmpStr = 0;
   signed char errSC = 0;
   unsigned long seqUL = 0; /*sequence on*/

   signed char pFullIdBl = def_pFullId_mainSamToAln;
   signed char *samFileStr = 0;
   signed char *refFileStr = 0;
   signed char *outFileStr = 0;

   struct alnSet alnSetStackST;

   struct seqST *refHeapAryST = 0;
   signed long lenRefSL = 0;
   signed long sizeRefSL = 0;
   signed long indexSL = 0;
 

   struct samEntry samStackST;

   FILE *samFILE = 0;
   FILE *outFILE = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^    - initialize, get user input, open files
   ^    o main sec02 sub01:
   ^      - initialize structures
   ^    o main sec02 sub02:
   ^      - get user input
   ^    o main sec02 sub03:
   ^      - set up structures
   ^    o main sec02 sub04:
   ^      - open sam file
   ^    o main sec02 sub05:
   ^      - open reference file and get reference
   ^    o main sec02 sub06:
   ^      - open output file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec02 Sub01:
   *   - initialize structures
   \*****************************************************/

   init_alnSet(&alnSetStackST);
   init_samEntry(&samStackST);

   alnSetStackST.lineWrapUI = def_lineWrap_mainSamToAln;
   alnSetStackST.pBasePosBl = def_ppos_mainSamToAln;
   alnSetStackST.pFullAlnBl = def_pmask_mainSamToAln;

   /*****************************************************\
   * Main Sec02 Sub02:
   *   - get user input
   \*****************************************************/

   errSC =
      input_mainSamToAln(
         numArgsSI,
         argAryStr,
         &samFileStr,
         &refFileStr,
         &outFileStr,
         &pFullIdBl,
         &alnSetStackST
      );

   if(errSC)
   { /*If: had an error*/
      --errSC;     /*convert help/version to 0*/
      errSC >>= 2; /*convert invalid input to 4*/
      goto cleanUp_main_sec04_sub04;
   } /*If: had an error*/

   /*****************************************************\
   * Main Sec02 Sub03:
   *   - set up structures
   \*****************************************************/

   setup_samEntry(&samStackST);

   /*****************************************************\
   * Main Sec02 Sub04:
   *   - open sam file
   \*****************************************************/

   if(! samFileStr || *samFileStr == '-')
      samFILE = stdin;

   else
   { /*Else: user provided a file*/
      samFILE =
         fopen(
            (char *) samFileStr,
            "r"
         );

      if(! samFILE)
      { /*If: could not open the sam file*/
         fprintf(
            stderr,
            "could not open -sam %s%s",
            samFileStr,
            str_endLine
         );

         goto fileErr_main_sec04_sub03;
      } /*If: could not open the sam file*/
   } /*Else: user provided a file*/

   /*****************************************************\
   * Main Sec02 Sub05:
   *   - open reference file and get reference
   *   o main sec02 sub05 cat01:
   *     - open reference file
   *   o main sec02 sub05 cat02:
   *     - get all sequences from fasta file
   *   o main sec02 sub05 cat03:
   *     - remove name white space and sort by id
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub05 Cat01:
   +   - open reference file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   outFILE =
      fopen(
         (char *) refFileStr,
         "r"
      );

   if(! outFILE)
   { /*If: could not open reference file*/
      fprintf(
         stderr,
         "could not open -ref %s%s",
         refFileStr,
         str_endLine
      );

      goto fileErr_main_sec04_sub03;
   } /*If: could not open reference file*/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub05 Cat02:
   +   - get all sequences from fasta file
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   refHeapAryST =
      readFaFile_seqST(
         outFILE,
         &lenRefSL,
         &sizeRefSL,
         &errSC
      ); /*get reference sequences from file*/

   if(errSC)
   { /*If: had error*/
      if(errSC & def_memErr_seqST)
      { /*If: memory error*/
         fprintf(
            stderr,
            "memory error reading -ref %s%s",
            refFileStr,
            str_endLine
         );

         goto memErr_main_sec04_sub02;
      } /*If: memory error*/

      else
      { /*If: not a fasta file*/
         fprintf(
            stderr,
            "-ref %s is not a fasta file%s",
            refFileStr,
            str_endLine
         );

         goto fileErr_main_sec04_sub03;
      } /*If: not a fasta file*/
   } /*If: had error*/

   fclose(outFILE);
   outFILE = 0;

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Main Sec02 Sub05 Cat03:
   +   - remove name white space and sort by id
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   for(
      indexSL = 0;
      indexSL < lenRefSL;
      ++indexSL
   ){ /*Loop: remove white space from names*/
      tmpStr = refHeapAryST[indexSL].idStr;
      tmpStr += endWhite_ulCp(tmpStr);
      *tmpStr = '\0';
   } /*Loop: remove white space from names*/

   sort_seqST(
      refHeapAryST,
      lenRefSL
   );

   /*****************************************************\
   * Main Sec02 Sub06:
   *   - open output file
   \*****************************************************/

   if(! outFileStr || *outFileStr == '-')
      outFILE = stdout;

   else
   { /*Else: user provided an output file*/
      outFILE =
         fopen(
            (char *) outFileStr,
            "w"
         );

      if(! outFILE)
      { /*If: could not open output file*/
         fprintf(
            stderr,
            "could not open -out %s%s",
            outFileStr,
            str_endLine
         );

         goto fileErr_main_sec04_sub03;
      } /*If: could not open output file*/
   } /*Else: user provided an output file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^    - print out alignments
   ^    o main sec03 sub01:
   ^      - get first sam file line
   ^    o main sec03 sub02:
   ^      - print global header (program history)
   ^    o main sec03 sub03:
   ^     - add mainSamToAln to program header
   ^    o main sec03 sub04:
   ^     - print alingments
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec03 Sub01:
   *   - get first sam file line
   \*****************************************************/

   errSC = get_samEntry(&samStackST, samFILE);

   if(errSC)
   { /*If: had an error*/
      if(errSC == def_EOF_samEntry)
      { /*If: had empty file*/
          fprintf(
             stderr,
             "-sam %s is empty%s",
             samFileStr,
             str_endLine
           );

         goto fileErr_main_sec04_sub03;
      } /*If: had empty file*/

      fprintf(
         stderr,
         "memory error reading first line of -sam %s%s",
         samFileStr,
         str_endLine
       );

      goto memErr_main_sec04_sub02;
   } /*If: had an error*/

   /*****************************************************\
   * Main Sec03 Sub02:
   *   - print global header (program history)
   \*****************************************************/

   errSC = psamPg_samToAln(&samStackST, samFILE, outFILE);

   if(errSC)
   { /*If: had an error*/
      if(errSC == def_EOF_samEntry)
      { /*If: had empty file*/
          fprintf(
             stderr,
             "-sam %s had only header entries%s",
             samFileStr,
             str_endLine
           );

         goto fileErr_main_sec04_sub03;
      } /*If: had empty file*/

      fprintf(
         stderr,
         "memory error header/first read of -sam %s%s",
         samFileStr,
         str_endLine
       );

      goto memErr_main_sec04_sub02;
   } /*If: had an error*/

   /*****************************************************\
   * Main Sec03 Sub03:
   *   - add mainSamToAln to program header
   \*****************************************************/

   fprintf(
     (FILE *) outFILE,
     "#PG\tID:%s\tPN:%s\tVN:bioTools_%i-%02i-%02i\tCL:%s",
     "samToAln",
     "samToAln",
     def_year_bioTools,
     def_month_bioTools,
     def_day_bioTools,
     "samToAln"
   ); /*print command without user arguments*/

   for(
      seqUL = 1;
      seqUL < (unsigned long) numArgsSI;
      ++seqUL
   ){ /*Loop: print out command arguments*/
      fprintf(
         (FILE *) outFILE,
         " %s",
         argAryStr[seqUL]
      );
   } /*Loop: print out command arguments*/

   fprintf(
      (FILE*) outFILE,
      "%s%s",
      str_endLine,
      str_endLine
    );

   /*****************************************************\
   * Main Sec03 Sub04:
   *   - print alingments
   \*****************************************************/

   seqUL = 0;

   while(! errSC)
   { /*Loop: print alignments in sam file*/
      ++seqUL;

      if(*samStackST.seqStr == '*')
         goto nextRead_main_sec03_sub04;

      if(*samStackST.cigTypeStr == '*')
         goto nextRead_main_sec03_sub04;

      if(samStackST.flagUS & 4)
         goto nextRead_main_sec03_sub04;

      indexSL =
         search_seqST(
            refHeapAryST,
            samStackST.refIdStr,
            lenRefSL
         ); /*see if can find reference*/

      if(indexSL < 0)
      { /*If: could not find reference*/
         fprintf(
            stderr,
            "-ref %s missing %s, skipping read %s%s",
            refFileStr,
            samStackST.refIdStr,
            samStackST.qryIdStr,
            str_endLine
         ); /*warn about missing reference*/

         goto nextRead_main_sec03_sub04;
      } /*If: could not find reference*/

      if(pFullIdBl)
      { /*If: user wanted full reference id*/
         tmpStr = refHeapAryST[indexSL].idStr;
         tmpStr += endWhite_ulCp(tmpStr);

         if(
               tmpStr - refHeapAryST[indexSL].idStr
            <  refHeapAryST[indexSL].idLenSL
         ) *tmpStr = ' '; /*Id was prematurely shortend*/
      } /*If: user wanted full reference id*/

      errSC =
         paln_samToAln(
            &samStackST,
            &refHeapAryST[indexSL],
            &alnSetStackST,
            outFILE
         ); /*print alignment*/

      if(pFullIdBl)
         *tmpStr = '\0';

      nextRead_main_sec03_sub04:;

      errSC = get_samEntry(&samStackST, samFILE);
   } /*Loop: print alignments in sam file*/

   if(errSC & def_memErr_samEntry)
   { /*If: had an error*/
      fprintf(
         stderr,
         "memory error for read %lu of -sam %s%s",
         seqUL,
         samFileStr,
         str_endLine
       );

      goto memErr_main_sec04_sub02;
   } /*If: had an error*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^   - clean up and return
   ^   o main sec04 sub01:
   ^     - no error clean up
   ^   o main sec04 sub02:
   ^     - memory error clean up
   ^   o main sec04 sub03:
   ^     - file error clean up
   ^   o main sec04 sub04:
   ^     - general clean up (error and no error)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Main Sec04 Sub01:
   *   - no error clean up
   \*****************************************************/

   errSC = 0;
   goto cleanUp_main_sec04_sub04;

   /*****************************************************\
   * Main Sec04 Sub02:
   *   - memory error clean up
   \*****************************************************/

   memErr_main_sec04_sub02:;
   errSC = 1;
   goto cleanUp_main_sec04_sub04;

   /*****************************************************\
   * Main Sec04 Sub03:
   *   - file error clean up
   \*****************************************************/

   fileErr_main_sec04_sub03:;
   errSC = 2;
   goto cleanUp_main_sec04_sub04;

   /*****************************************************\
   * Main Sec04 Sub04:
   *   - general clean up (error and no error)
   \*****************************************************/

   cleanUp_main_sec04_sub04:;

   freeStack_alnSet(&alnSetStackST);
   freeStack_samEntry(&samStackST);


   if(refHeapAryST)
      freeHeapAry_seqST(
         refHeapAryST,
         sizeRefSL
      );
   refHeapAryST = 0;

   if(
         samFILE
      && samFILE != stdin
      && samFILE != stdout
   ) fclose(samFILE);
   samFILE = 0;


   if(
         outFILE
      && outFILE != stdin
      && outFILE != stdout
   ) fclose(outFILE);
   outFILE = 0;


   return(errSC);
} /*main*/

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However,for cases were the public domain is not
:   suitable,such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible,this code is under the
:   MIT license.
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy,modify,publish,use,compile,
:   sell,or distribute this software,either in source
:   code form or as a compiled binary,for any purpose,
:   commercial or non-commercial,and by any means.
: 
: In jurisdictions that recognize copyright laws,the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND,EXPRESS OR IMPLIED,INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY,WHETHER IN AN ACTION OF
:   CONTRACT,TORT OR OTHERWISE,ARISING FROM,OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information,please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2024 jeremyButtler
: 
: Permission is hereby granted,free of charge,to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"),to
:   deal in the Software without restriction,including
:   without limitation the rights to use,copy,modify,
:   merge,publish,distribute,sublicense,and/or sell
:   copies of the Software,and to permit persons to whom
:   the Software is furnished to do so,subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND,EXPRESS OR IMPLIED,INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY,FITNESS
:   FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM,DAMAGES OR OTHER LIABILITY,WHETHER IN
:   AN ACTION OF CONTRACT,TORT OR OTHERWISE,ARISING
:   FROM,OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
