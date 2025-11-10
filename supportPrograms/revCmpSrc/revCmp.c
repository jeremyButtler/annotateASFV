/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' revCmp SOF: Start Of File
'   - reverse complements all sequences in a fasta file
'   o header:
'     - included libraries
'   o fun01: pversion_revCmp
'     - Prints version number
'   o fun02 pHelpMsg:
'     - Prints out the help message
'   o main main:
'     - Reverse complements all sequences in a fasta file
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

#include "../genLib/ulCp.h"
#include "../genBio/seqST.h"

/*.h files only*/
#include "../genLib/endLine.h"
#include "../bioTools.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - .c  #include "../genLib/fileFun.h"
!   - .c  #include "../genLib/ulCp.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*--------------------------------------------------------\
| Fun01: pversion_revCmp
|  - Prints version number
| Input:
|  - outFILE:
|    o Pointer to FILE to print to
| Output:
|  - Prints:
|    o version number to outFILE
\--------------------------------------------------------*/
void
pversion_revCmp(
   void *outFILE
){
   fprintf(
     (FILE *) outFILE,
     "revCmp from bioTools version: %i-%02i-%02i%s",
     def_year_bioTools,
     def_month_bioTools,
     def_day_bioTools,
     str_endLine
   );
} /*pversion_revCmp*/


/*--------------------------------------------------------\
| Fun02: phelp_revCmp
|  - Prints out the help message
| Input:
|  - outFILE:
|    o Pointer to FILE to print to
| Output:
|  - Prints:
|    o help message to outFILE
\--------------------------------------------------------*/
void
phelp_revCmp(
   void *outFILE
){
   fprintf(
       (FILE *) outFILE,
       "revCmpSeq -out out.fa file_1.fa file_2.fa ...%s",
       str_endLine
   );

   fprintf(
     (FILE *) outFILE,
     "  - reverse complement sequences in fasta a file%s",
     str_endLine
   );


   fprintf((FILE *) outFILE, "Input:%s", str_endLine);
   fprintf(
     (FILE *) outFILE,
     "  file.fa: [Required]%s",
     str_endLine
   );

   fprintf(
     (FILE *) outFILE,
     "    o Fasta file with sequences to reverse%s",
     str_endLine
   );

   fprintf(
     (FILE *) outFILE,
     "      complement%s",
     str_endLine
   );


   fprintf(
     (FILE *) outFILE,
     "  -out out.fa: [Optional: stdout]%s",
     str_endLine
   );
   fprintf(
     (FILE *) outFILE,
     "    o File to save reverse complement seqs to%s",
     str_endLine
   );
   fprintf(
     (FILE *) outFILE,
     "    o Use `-` for stdout%s",
     str_endLine
   );
   fprintf(
     (FILE *) outFILE,
     "    o This must be the first arugment%s",
     str_endLine
   );


   fprintf((FILE *) outFILE, "Output:%s", str_endLine);
   fprintf(
     (FILE *) outFILE,
     "  - Reverse complement sequences to stdout%s",
     str_endLine
   );
} /*phelp_revCmp*/

/*--------------------------------------------------------\
| Main: main
|  - Reverse complements the input sequence
| Input:
|  - numArgsSI:
|    o Number of arguments the user input
|  - argAryStr:
|    o Arguments and parameters the user input
| Output:
|  - Prints:
|    o reverse complement sequence to stdout
\--------------------------------------------------------*/
int main(
   int numArgsSI,   /*Number of arguments the user input*/
   char *argAryStr[] /*Arguments/parametersminput by user*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Main: main
   '  - Reverse complements sequences in argument 1
   '  o main sec01:
   '    - Variable declerations
   '  o main sec02:
   '    - Check user input and initialize variables
   '  o main sec03:
   '    - Get and reverse complement fasta entries
   '  o main sec04:
   '    - Clean up and exit
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec01:
   ^  - Variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;
   signed int numSeqSI = 0;
   signed int siArg = 1;
   signed char stdinBl = 0;

   FILE *inFILE = 0;
   FILE *outFILE = stdout;
   struct seqST seqStackST;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec02:
   ^  - Check user input and initialize variables
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   init_seqST(&seqStackST);

   if(numArgsSI < 2)
   { /*If: printing help message*/
      phelp_revCmp(stdout);
      goto clean_main_sec04;
   } /*If: printing help message*/


   if(argAryStr[1][0] == '-')
   { /*If: help message, output file or version number*/
      if(argAryStr[1][1] == '-')
      { /*If: used two dashes*/
         if(argAryStr[1][2] == 'h')
         { /*Else If: likely help message*/
            phelp_revCmp(stdout);
            goto clean_main_sec04;
         } /*Else If: likely help message*/

         else if(argAryStr[1][2] == 'v')
         { /*Else If: likely version number*/
            pversion_revCmp(stdout);
            goto clean_main_sec04;
         } /*Else If: likely version number*/

         else if(argAryStr[1][2] == 'o')
         { /*Else If: likley output file*/
            if(argAryStr[2][0] != '-')
            { /*If: file was input*/
               outFILE = 0;
               outFILE = fopen(argAryStr[2], "w");
               if(! outFILE)
                  goto outFileErr_main_sec04;
            } /*If: file was input*/
 
            siArg = 3;
         } /*Else If: likely output file*/

         else
            goto unkownInput_main_sec04; /*no idea*/
      } /*If: used two dashes*/

      else if(! argAryStr[1][1])
         ; /*uer wants stdin input*/

      else if(argAryStr[1][1] == 'h')
      { /*Else If: user likely requested help message*/
         phelp_revCmp(stdout);
         goto clean_main_sec04;
      } /*Else If: user likely requested help message*/

      else if(argAryStr[1][1] == 'v')
      { /*Else If: user likely requested version number*/
         pversion_revCmp(stdout);
         goto clean_main_sec04;
      } /*Else If: user likely requested version number*/

      else if(argAryStr[1][1] == 'o')
      { /*Else If: user input output file*/
         if(argAryStr[2][0] != '-')
         { /*If: file was input*/
            outFILE = 0;
            outFILE = fopen(argAryStr[2], "w");
            if(! outFILE)
               goto outFileErr_main_sec04;
         } /*If: file was input*/
 
         siArg = 3;
      } /*Else If: user input output file*/

      else
         goto unkownInput_main_sec04; /*no idea what is*/
   } /*If: help message, output file or version number*/


   if(
      ! eqlNull_ulCp(
         (signed char *) "help",
         (signed char *) argAryStr[1]
      )
   ){ /*If: the help message was requested*/
      phelp_revCmp(stdout);
      goto clean_main_sec04;
   } /*If: the help message was requested*/

   else if(
      ! eqlNull_ulCp(
         (signed char *) "version",
         (signed char *) argAryStr[1]
      )
   ){ /*Else If: version number requested*/
      pversion_revCmp(stdout);
      goto clean_main_sec04;
   } /*Else If: version number requested*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec03:
   ^  - Get and reverse complement fasta entries
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(siArg < numArgsSI)
   { /*Loop: got through fasta files*/
      if(! inFILE) ;
      else if(inFILE == stdin) ;
      else
         fclose(inFILE);
      inFILE = 0;

      if(argAryStr[siArg][0] == '-')
      { /*If: stdin input*/
         if(argAryStr[siArg][1])
            goto unkownInput_main_sec04;

         if(stdinBl)
         { /*If: mulitple stdin inputs requested*/
            ++siArg;
            continue;
         } /*If: mulitple stdin inputs requested*/

         inFILE = stdin;
         stdinBl = 1;
      } /*If: stdin input*/

      else
      { /*Else: need to open the file*/
         inFILE = fopen(argAryStr[siArg], "r");
         if(! inFILE)
            goto faOpenErr_main_sec04;
      } /*Else: need to open the file*/

      blank_seqST(&seqStackST);
      errSC = getFa_seqST(inFILE, &seqStackST);
      numSeqSI = 0;

      if(errSC == def_memErr_seqST)
         goto memErr_main_sec04;
      else if(errSC)
         goto faFileErr_main_sec04;

      while(! errSC)
      { /*Loop: though all sequences in the fasta file*/
         ++numSeqSI;

         revComp_seqST(&seqStackST);

         fprintf(
            outFILE,
            ">%s%s%s%s",
            seqStackST.idStr,
            str_endLine,
            seqStackST.seqStr,
            str_endLine
         ); /*Print out the sequence*/

         blank_seqST(&seqStackST);
         errSC = getFa_seqST(inFILE, &seqStackST);
      } /*LoopCheck user input*/

      if(errSC == def_memErr_seqST)
         goto memErr_main_sec04;
      else if(errSC != def_EOF_seqST)
         goto faFileErr_main_sec04;
      else
         errSC = 0;

      ++siArg;
   } /*Loop: got through fasta files*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Main Sec04:
   ^  - Clean up and exit
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto clean_main_sec04;

   memErr_main_sec04:;
      fprintf(stderr, "memory error%s", str_endLine);
      errSC = 1;
      goto clean_main_sec04;

   outFileErr_main_sec04:;
      fprintf(
         stderr,
         "failed to open -out %s%s",
         argAryStr[2], /*-out is 1st value, file is 3rd*/
         str_endLine

      );
      errSC = 2;
      goto clean_main_sec04;

   faOpenErr_main_sec04:;
      fprintf(
         stderr,
         "failed to open %s%s",
         argAryStr[siArg],
         str_endLine

      );
      errSC = 3;
      goto clean_main_sec04;

   faFileErr_main_sec04:;
      fprintf(
         stderr,
         "line %i in %s is not a fasta entry%s",
         numSeqSI,
         argAryStr[siArg],
         str_endLine

      );
      errSC = 4;
      goto clean_main_sec04;

   unkownInput_main_sec04:;
      fprintf(
         stderr,
         "%s is not -out out.fa%s",
         argAryStr[siArg],
         str_endLine

      );
      errSC = 4;
      goto clean_main_sec04;

   clean_main_sec04:;
      freeStack_seqST(&seqStackST);

      if(! inFILE) ;
      else if(inFILE == stdin) ;
      else if(inFILE == stdout) ;
      else if(inFILE == stderr) ;
      else fclose(inFILE);
      inFILE = 0;

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
