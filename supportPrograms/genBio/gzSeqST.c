/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' seqST SOF: Start Of File
'   - holds gz functions for seqST
'   - the reason this is separate is because my inflate
'     libraries adds around 18kb to all programs, even
'     ones that do not use it
'   o header:
'     - included libraries
'   o fun01: get_gzSeqST
'     - gets a sequence from a fastx/fastx.gz file
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>
#include "gzSeqST.h"
#include "seqST.h"
#include "../genLib/ulCp.h"
#include "../genLib/fileFun.h"
#include "../genLib/inflate.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - .c  #include "../genLib/checkSum.h"
!   - .c  #include "../genLib/genMath.h"
!   - .c  #include "../genLib/endin.h"
!   - .h  #include "../genLib/endLine.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun01: get_gzSeqST
|   - gets a sequence from a fastx/fastx.gz file
| Input:
|   - fileSTPtr:
|     o file_inflate struct to use in uncompression
|   - typeSCPtr:
|     o is set to or has the correct file type
|   - seqSTPtr:
|     o seqST struct pointer to add new sequence to
|   - inFILE:
|     o FILE pointer to file to add (new file)
|     o 0 if already called get_gzSeqST for this file
|     o must be opened in binary format
|       * fopen(<name>, "rb");
|     o you should set inFILE to 0/null after calling,
|       because freeStack_file_inflate 
|       freeHeap_file_inflate will close inFILE for you
| Output:
|   - Modifies:
|     o if have inFILE, sets fileSTPtr be blanked and then
|       have inFILE
|     o if have inFILE, typeSCPtr to have file type 
|       o def_noType_gzSeqST if end of file
|       o def_faType_gzSeqST if fasta file
|       o def_faType_gzSeqST | def_gzType_gzSeqST if
|         fasta.gz file
|       o def_fqType_gzSeqST if fastq file
|       o def_fqType_gzSeqST | def_gzType_gzSeqST if
|         fastq.gz file
|     o fileSTPtr to be on next read/sequence
|     o seqSTPtr to have next read/sequence
|   - Returns:
|     o 0 for no errors
|     o def_EOF_seqST if at end of file
|     o def_memErr_seqST for memory errors
|     o def_badLine_seqST | def_fileErr_seqST for invalid
|       file entries
|     o def_fileErr_seqST for file errors (including gzip
|       checks)
\-------------------------------------------------------*/
signed char
get_gzSeqST(
   struct file_inflate *fileSTPtr, /*has file to extract*/
   signed char *typeSCPtr,         /*gets/has file type*/
   struct seqST *seqSTPtr,         /*gets sequence*/
   void *inFILE
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun01 TOC:
   '   - gets a sequence from a fastx/fastx.gz file
   '   o fun01 sec01:
   '     - variable declarations
   '   o fun01 sec02:
   '     - add input file (if one) to file_inflate struct
   '   o fun01 sec03:
   '     - get read id
   '   o fun01 sec04:
   '     - get sequence
   '   o fun01 sec05:
   '     - get q-score entry (fastq only)
   '   o fun01 sec06:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #define extraBytes_fun01 1024;
   signed char errSC = 0;
   signed char *tmpStr = 0;
   signed char tmpSC = 0;
   signed char tmpArySC[4];

   signed long lenSL = 0;
   signed long bytesSL = 0;
   signed long tmpSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec02:
   ^   - add input file (if one) to file_inflate struct
   ^   o fun01 sec02 sub01:
   ^     - get memory for id (if no memory)
   ^   o fun01 sec02 sub02:
   ^     - if file input; get file type and id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec02 Sub01:
   *   - get memory for id (if no memory)
   \*****************************************************/

   if(! seqSTPtr->idStr)
   { /*If: need memory for read id*/
      seqSTPtr->idStr = malloc(136 * sizeof(signed char));
      if(! seqSTPtr)
         goto memErr_fun01_sec06;
      seqSTPtr->idSizeSL = 128;
   } /*If: need memory for read id*/

   blank_seqST(seqSTPtr);

   /*****************************************************\
   * Fun01 Sec02 Sub02:
   *   - if file input; get file type and id
   *   o fun01 sec02 sub02 cat01:
   *     - find if .gz file and get first line
   *   o fun01 sec02 sub02 cat02:
   *     - check if .fastq file .fasta
   *   o fun01 sec02 sub02 cat03:
   *     - first line empty, scan for fasta header
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun01 Sec02 Sub02 Cat01:
   +   - find if .gz file and get first line
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(inFILE)
   { /*If: user is starting new file*/
      *typeSCPtr = 0;

      if(
        ! fread(
           (char *) seqSTPtr->idStr,
           sizeof(signed char),
           2,
           inFILE
         )
      ) goto fileErr_fun01_sec06;

      if(seqSTPtr->idStr[0] == 0x1f) /*0x1f = 31*/
      { /*If: gz file or invalid file*/
         *typeSCPtr = def_gzType_gzSeqST;

         errSC =
           add_file_inflate(
              inFILE,
              fileSTPtr,
              seqSTPtr->idStr
           );

         if(errSC == def_memErr_inflate)
            goto memErr_fun01_sec06;
         else if(errSC)
            goto fileErr_fun01_sec06;

         seqSTPtr->idLenSL =
            get_inflate(
               fileSTPtr,
               seqSTPtr->idSizeSL,
               1, /*till line ends or buffer is full*/
               seqSTPtr->idStr,
               &errSC
            ); /*get first line or first 128 characters*/

         if(errSC == def_memErr_inflate)
            goto memErr_fun01_sec06;
         else if(errSC)
            goto fileErr_fun01_sec06;
      } /*If: gz file or invalid file*/

      else
      { /*Else: file is not gz compressed*/
         blank_file_inflate(fileSTPtr, 2);
         fileSTPtr->zipFILE = inFILE;

         seqSTPtr->idLenSL =
            getLine_fileFun(
               fileSTPtr->zipFILE,
               &seqSTPtr->idStr[2],
               seqSTPtr->idSizeSL - 2,
               &tmpSL
            ); /*get first line into buffer*/

         if(! seqSTPtr->idLenSL)
            goto fileErr_fun01_sec06;
         seqSTPtr->idLenSL += 2;
      } /*Else: file is not gz compressed*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun01 Sec02 Sub02 Cat02:
      +   - check if .fastq file .fasta
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(seqSTPtr->idStr[0] == '@')
         *typeSCPtr |= def_fqType_gzSeqST;
      else if(seqSTPtr->idStr[0] == '>')
         *typeSCPtr |= def_faType_gzSeqST;
      else if(seqSTPtr->idStr[1] == '>')
         *typeSCPtr |= def_faType_gzSeqST;
      else if(seqSTPtr->idStr[2] == '>')
         *typeSCPtr |= def_faType_gzSeqST;
      else
      { /*Else: may need to get past blank lines*/
         if(seqSTPtr->idStr[0] > 32)
            goto fileErr_fun01_sec06; /*unkown file*/
         else if(seqSTPtr->idStr[1] > 32)
            goto fileErr_fun01_sec06; /*unkown file*/
         else if(seqSTPtr->idStr[2] > 32)
            goto fileErr_fun01_sec06; /*unkown file*/

         if(seqSTPtr->idLenSL > 3)
         { /*If: have something to copy*/
            seqSTPtr->idLenSL -= 3;
            cpLen_ulCp(
               seqSTPtr->idStr,
               &seqSTPtr->idStr[3],
               seqSTPtr->idLenSL
            );
         } /*If: have something to copy*/

         else
            seqSTPtr->idLenSL = 0;

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun01 Sec02 Sub02 Cat03:
         +   - first line empty, scan for fasta header
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         while(! (*typeSCPtr & def_faType_gzSeqST) )
         { /*Loop: check if fasta or invalid file*/
            tmpStr = seqSTPtr->idStr;
            while(*tmpStr < 33 && *tmpStr)
               ++tmpStr; /*find first non-white space*/
                
            if(*tmpStr == '>')
            { /*If: found file type*/
               *typeSCPtr |= def_faType_gzSeqST;
               break;
            } /*If: found file type*/

            else if(*tmpStr)
               goto fileErr_fun01_sec06; /*not .fasta*/
            else
               seqSTPtr->idLenSL = 0; /*only white space*/


            if(*typeSCPtr & def_gzType_gzSeqST)
            { /*If: reading gz file*/
               seqSTPtr->idLenSL =
                  get_inflate(
                     fileSTPtr,
                     seqSTPtr->idSizeSL,
                     1, /*one line or full buffer*/
                     seqSTPtr->idStr,
                     &errSC
                  );

               if(errSC == def_memErr_inflate)
                  goto memErr_fun01_sec06;
               else if(errSC)
                  goto fileErr_fun01_sec06;
            } /*If: reading gz file*/


            else
            { /*Else: fasta file (uncompressed)*/
               seqSTPtr->idLenSL =
                  getLine_fileFun(
                     fileSTPtr->zipFILE,
                     seqSTPtr->idStr,
                     seqSTPtr->idSizeSL,
                     &tmpSL
                  ); /*get first line into buffer*/

               if(! seqSTPtr->idLenSL)
                  goto fileErr_fun01_sec06;
            } /*Else: fasta file (uncompressed)*/
         } /*Loop: check if fasta or invalid file*/
      } /*Else: may need to get past blank lines*/
   } /*If: user is starting new file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec03:
   ^   - get read id
   ^   o fun01 sec03 sub01:
   ^     - get first line of read id (if no input file)
   ^   o fun01 sec03 sub02:
   ^     - get read id
   ^   o fun01 sec03 sub03:
   ^     - remove line break and header marker from id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec03 Sub01:
   *   - get first line of read id (if no input file)
   \*****************************************************/

   if(! inFILE)
   { /*If: need to get the first line still*/

      if(*typeSCPtr & def_gzType_gzSeqST)
      { /*If: is gz file*/
         errSC = 0;

         while(! errSC)
         { /*Loop: get next line*/
            seqSTPtr->idLenSL =
               get_inflate(
                  fileSTPtr,
                  seqSTPtr->idSizeSL,
                  1, /*till line ends or buffer is full*/
                  seqSTPtr->idStr,
                  &errSC
               ); /*get first line*/

            if(errSC == def_memErr_inflate)
               goto memErr_fun01_sec06;

            else if(errSC == def_eof_inflate)
            { /*Else If: at end of gz file or entry*/
               errSC = next_file_inflate(fileSTPtr, 0, 0);

               if(errSC == def_eof_inflate)
                  goto eof_fun01_sec06;
               else if(errSC)
                  goto fileErr_fun01_sec06;
            } /*Else If: at end of gz file or entry*/

            else
               break; /*got one line*/
         } /*Loop: get next line*/
      } /*If: is gz file*/

      else
      { /*Else: reading uncompressed file*/
         seqSTPtr->idLenSL =
            getLine_fileFun(
               fileSTPtr->zipFILE,
               seqSTPtr->idStr,
               seqSTPtr->idSizeSL,
               &tmpSL
            ); /*get first line into buffer*/
        if(! seqSTPtr->idLenSL)
           goto eof_fun01_sec06;
      } /*Else: reading uncompressed file*/
   } /*If: need to get the first line still*/

   /*****************************************************\
   * Fun01 Sec03 Sub02:
   *   - get read id
   \*****************************************************/

   while(
         seqSTPtr->idStr[seqSTPtr->idLenSL - 1] != '\n'
      && seqSTPtr->idStr[seqSTPtr->idLenSL - 1] != '\r'
   ){ /*Loop: get read id*/
      tmpStr =
         realloc(
            seqSTPtr->idStr,
            (seqSTPtr->idLenSL + 136)
               * sizeof(signed char)
         );

      if(! tmpStr)
         goto memErr_fun01_sec06;
      seqSTPtr->idStr = tmpStr;
      seqSTPtr->idSizeSL += 128;

      if(*typeSCPtr & def_gzType_gzSeqST)
      { /*If: is gz file*/
         seqSTPtr->idLenSL +=
            get_inflate(
               fileSTPtr,
               128,
               1, /*till line ends or buffer is full*/
               &seqSTPtr->idStr[seqSTPtr->idLenSL],
               &errSC
            );

            if(errSC == def_memErr_inflate)
               goto memErr_fun01_sec06;
            else if(errSC == def_eof_inflate)
               goto badLine_fun01_sec06;
            else if(errSC)
               goto fileErr_fun01_sec06;
      } /*If: is gz file*/

      else
      { /*Else: uncompressed file*/
         seqSTPtr->idLenSL +=
            getLine_fileFun(
               fileSTPtr->zipFILE,
               &seqSTPtr->idStr[seqSTPtr->idLenSL],
               128,
               &tmpSL
            ); /*get first line into buffer*/
        if(! seqSTPtr->idLenSL)
           goto badLine_fun01_sec06;
      } /*Else: uncompressed file*/
   } /*Loop: get read id*/

   /*****************************************************\
   * Fun01 Sec03 Sub03:
   *   - remove line break and header marker from id
   \*****************************************************/

   while(
         seqSTPtr->idStr[seqSTPtr->idLenSL - 1] == '\n'
      || seqSTPtr->idStr[seqSTPtr->idLenSL - 1] == '\r'
   ) --seqSTPtr->idLenSL;

   tmpSL = 0;
   while(
         seqSTPtr->idStr[tmpSL] < 33
      && tmpSL < seqSTPtr->idLenSL
   ) ++tmpSL;

   if(seqSTPtr->idStr[tmpSL] == '@')
      ++tmpSL;
   else if(seqSTPtr->idStr[tmpSL] == '>')
      ++tmpSL;

   cpLen_ulCp(
      seqSTPtr->idStr,
      &seqSTPtr->idStr[tmpSL], /*move past header symbol*/
      seqSTPtr->idLenSL - tmpSL
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec04:
   ^   - get sequence
   ^   o fun01 sec04 sub01:
   ^     - get first line or part of line in sequence
   ^   o fun01 sec04 sub02:
   ^     - check if reached end of fastq sequence
   ^   o fun01 sec04 sub03:
   ^     - check if reached end of fasta sequence
   ^   o fun01 sec04 sub04:
   ^     - make sure have memory for next line
   ^   o fun01 sec04 sub05:
   ^     - read in next line or part of line
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec04 Sub01:
   *   - get first line or part of line in sequence
   \*****************************************************/

   if(! seqSTPtr->seqStr)
   { /*If: need to get memory for sequence*/
      seqSTPtr->seqStr =
         malloc(1032 * sizeof(signed char));
      if(! seqSTPtr->seqStr)
         goto memErr_fun01_sec06;
      seqSTPtr->seqSizeSL = 1024;
   } /*If: need to get memory for sequence*/


   if(*typeSCPtr & def_gzType_gzSeqST)
   { /*If: reading gz file*/
      lenSL =
         get_inflate(
            fileSTPtr,
            seqSTPtr->seqSizeSL,
            1, /*till line ends or buffer is full*/
            seqSTPtr->seqStr,
            &errSC
         ); /*get first sequence line*/

      if(errSC == def_memErr_inflate)
         goto memErr_fun01_sec06;
      else if(errSC == def_eof_inflate)
         goto badLine_fun01_sec06;
      else if(errSC)
         goto fileErr_fun01_sec06;
   } /*If: reading gz file*/

   else
   { /*Else: reading uncompressed file*/
      lenSL =
         getLine_fileFun(
            fileSTPtr->zipFILE,
            seqSTPtr->seqStr,
            seqSTPtr->seqSizeSL,
            &tmpSL
         );

      if(! lenSL)
         goto badLine_fun01_sec06;
   } /*Else: reading uncompressed file*/

   /*****************************************************\
   * Fun01 Sec04 Sub02:
   *   - check if reached end of fastq sequence
   \*****************************************************/

   while(1 == 1)
   { /*Loop: read in sequence*/
      tmpSC = 0;
      tmpArySC[0] = 0;

      if(*typeSCPtr & def_fqType_gzSeqST)
      { /*If: looking for fastq file ending*/
         if(seqSTPtr->seqStr[seqSTPtr->seqLenSL] == '+')
            break;
         seqSTPtr->seqLenSL += lenSL;

         if(*typeSCPtr & def_gzType_gzSeqST)
            bytesSL += lenSL;
         else
            bytesSL += tmpSL;
      } /*If: looking for fastq file ending*/

      /**************************************************\
      * Fun01 Sec04 Sub03:
      *   - check if reached end of fasta sequence
      *   o  fun01 sec04 sub03 cat01:
      *     - next fasta char from .gz file + .fa* check
      *   o fun01 sec04 sub03 cat02:
      *     - get next fasta char from .fa file
      *   o fun01 sec04 sub03 cat03:
      *     - check if have header or sequence
      \**************************************************/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun01 Sec04 Sub03 Cat01:
      +   - get next fasta char from .gz file + .fa* check
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else if(*typeSCPtr & def_faType_gzSeqST)
      { /*Else If: fasta file*/
         seqSTPtr->seqLenSL += lenSL;

         if(*typeSCPtr & def_gzType_gzSeqST)
         { /*If: reading from gz file*/
            while(! errSC)
            { /*Loop: find if next line is header*/
               lenSL =
                  getc_file_inflate(fileSTPtr, tmpArySC);

               if(lenSL < 0)
               { /*If: error*/
                  lenSL *= -1;

                  if(lenSL == def_memErr_inflate)
                     goto memErr_fun01_sec06;
                  else if(lenSL == def_eof_inflate)
                     goto done_fun01_sec06;
                     /*call eof on next sequence read*/
                  else
                     goto fileErr_fun01_sec06;
               } /*If: error*/

               tmpSC = tmpArySC[0];

               if(tmpSC > 32)
                  break;
            } /*Loop: find if next line is header*/
         } /*If: reading from gz file*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun01 Sec04 Sub03 Cat02:
         +   - get next fasta char from .fa file
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         else
         { /*Else: reading uncompressed file*/
            tmpSC  = 0;

            while(tmpSC != EOF)
            { /*Loop: find first non-white space*/
               tmpSC = fgetc(fileSTPtr->zipFILE);

               if(tmpSC > 32)
                  break;
               else if(tmpSC == EOF)
                  goto done_fun01_sec06;
            }  /*Loop: find first non-white space*/
         } /*Else: reading uncompressed file*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun01 Sec04 Sub03 Cat03:
         +   - check if have header or sequence
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         if(tmpSC == '>')
         { /*If: on next header*/
            if(*typeSCPtr & def_gzType_gzSeqST)
               unget_file_inflate(fileSTPtr, tmpSC);
            else
               ungetc(tmpSC, fileSTPtr->zipFILE);
            goto done_fun01_sec06;
         } /*If: on next header*/

         else if( (tmpSC & ~32) > 90 )
            goto fileErr_fun01_sec06; /*non a-z*/
         else if( (tmpSC & ~32) > 64 )
            ; /*a-z*/
         else if(tmpSC == '-')
            ; /*gap*/
         else
            goto badLine_fun01_sec06;
            /*not part of sequence*/
      } /*Else If: fasta file*/

      /**************************************************\
      * Fun01 Sec04 Sub04:
      *   - make sure have memory for next line
      \**************************************************/

      if(
            seqSTPtr->seqLenSL
         >= seqSTPtr->seqSizeSL - 128
      ){ /*If: need more memory*/
         seqSTPtr->seqSizeSL +=
            (seqSTPtr->seqSizeSL >> 1);

         tmpStr =
            realloc(
               seqSTPtr->seqStr,
               (seqSTPtr->seqSizeSL + 8)
                  * sizeof(signed char)
            );
         if(! tmpStr)
            goto memErr_fun01_sec06;
         seqSTPtr->seqStr = tmpStr;
      } /*If: need more memory*/

      /**************************************************\
      * Fun01 Sec04 Sub05:
      *   - read in next line or part of line
      \**************************************************/

      if(tmpSC)
         seqSTPtr->seqStr[seqSTPtr->seqLenSL++] = tmpSC;


      if(*typeSCPtr & def_gzType_gzSeqST)
      { /*If: reading gz file*/
         lenSL =
            get_inflate(
               fileSTPtr,
               seqSTPtr->seqSizeSL - seqSTPtr->seqLenSL,
               1, /*till line ends or buffer is full*/
               &seqSTPtr->seqStr[seqSTPtr->seqLenSL],
               &errSC
            );
         tmpSL = lenSL;

         if(errSC == def_memErr_inflate)
            goto memErr_fun01_sec06;
         else if(errSC == def_eof_inflate)
            goto eof_fun01_sec06;
         else if(errSC)
            goto fileErr_fun01_sec06;
      } /*If: reading gz file*/

      else
      { /*Else: reading uncompressed file*/
         lenSL =
            getLine_fileFun(
               fileSTPtr->zipFILE,
               &seqSTPtr->seqStr[seqSTPtr->seqLenSL],
               seqSTPtr->seqSizeSL - seqSTPtr->seqLenSL,
               &tmpSL
            );
            
         if(! lenSL && *typeSCPtr & def_fqType_gzSeqST)
            goto fileErr_fun01_sec06;
         else if(! lenSL)
            goto eof_fun01_sec06;
      } /*Else: reading uncompressed file*/
   } /*Loop: read in sequence*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec05:
   ^   - get q-score entry (fastq only)
   ^   o fun01 sec05 sub01:
   ^     - find length of q-score entry from sequence
   ^   o fun01 sec05 sub02:
   ^     - make sure have fastq spacer line read in
   ^   o fun01 sec05 sub03:
   ^     - get q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun01 Sec05 Sub01:
   *   - find length of q-score entry from sequence entry
   \*****************************************************/

   seqSTPtr->seqStr[seqSTPtr->seqLenSL] = 0;

   if(! seqSTPtr->qStr) ;

   else if(seqSTPtr->qSizeSL < bytesSL + 1)
   { /*Else If: q-score buffer is to small*/
      free(seqSTPtr->qStr);
      seqSTPtr->qStr = 0;
   } /*Else If: q-score buffer is to small*/

   if(! seqSTPtr->qStr)
   { /*If: need to get memory for sequence*/
      seqSTPtr->qStr =
         malloc((bytesSL + 8) * sizeof(signed char));
      if(! seqSTPtr->qStr)
         goto memErr_fun01_sec06;
      seqSTPtr->qSizeSL = bytesSL;
   } /*If: need to get memory for sequence*/

   /*****************************************************\
   * Fun01 Sec05 Sub02:
   *   - make sure have fastq spacer line read in
   \*****************************************************/

   if(
         seqSTPtr->seqStr[seqSTPtr->seqLenSL + lenSL - 1]
      == '\n'
   ) ;

   else if(
         seqSTPtr->seqStr[seqSTPtr->seqLenSL + lenSL - 1]
      == '\r'
   ) ;

   else
   { /*Else: need to finish reading in + entry*/
      do{ /*Loop: get full + line*/
         if(*typeSCPtr & def_gzType_gzSeqST)
         { /*If: reading gz file*/
            lenSL =
               get_inflate(
                  fileSTPtr,
                  (signed long) seqSTPtr->qSizeSL,
                  1, /*get line + have extra for break*/
                  seqSTPtr->qStr,
                  &errSC
               ); /*get first line*/

            if(errSC == def_memErr_inflate)
               goto memErr_fun01_sec06;
            if(errSC == def_eof_inflate)
               goto badLine_fun01_sec06;
            else if(errSC)
               goto fileErr_fun01_sec06;
         } /*If: reading gz file*/

         else
         { /*Else: uncompressed file*/
            lenSL =
               getLine_fileFun(
                  fileSTPtr->zipFILE,
                  seqSTPtr->qStr,
                  seqSTPtr->qSizeSL,
                  &tmpSL
            );
            if(! lenSL)
               goto badLine_fun01_sec06;
         } /*Else: uncompressed file*/

         tmpStr = seqSTPtr->qStr + lenSL - 1;
      } while(*tmpStr != '\n' && *tmpStr != '\0');
        /*Loop: get full + line*/
   } /*Else: need to finish reading in + entry*/

   /*****************************************************\
   * Fun01 Sec05 Sub03:
   *   - get q-score entry
   \*****************************************************/

   if(*typeSCPtr & def_gzType_gzSeqST)
   { /*If: reading gz file*/
      seqSTPtr->qLenSL +=
         get_inflate(
            fileSTPtr,
            bytesSL,
            4, /*buffer has extra bytes for line ending*/
            seqSTPtr->qStr,
            &errSC
         ); /*get first sequence line*/

      if(errSC == def_memErr_inflate)
         goto memErr_fun01_sec06;
      else if(errSC == def_eof_inflate)
         goto eof_fun01_sec06;
      else if(errSC)
         goto fileErr_fun01_sec06;
   } /*If: reading gz file*/

   else
   { /*Else: reading text fastq file*/
      seqSTPtr->qLenSL =
         fread(
            (char *) seqSTPtr->qStr,
            sizeof(signed char),
            bytesSL,
            fileSTPtr->zipFILE
         );
   } /*Else: reading text fastq file*/

   seqSTPtr->qStr[seqSTPtr->seqLenSL] = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun01 Sec06:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   done_fun01_sec06:;
      errSC = 0;
      goto noErr_fun01_sec06;

   eof_fun01_sec06:;
      errSC = def_EOF_seqST;
      if(*typeSCPtr & def_gzType_gzSeqST)
      { /*If: gz file, might have more entries*/
         errSC = next_file_inflate(fileSTPtr, 0, 0);
         if(! errSC)
            goto done_fun01_sec06;
         else if(errSC != def_eof_inflate)
            goto fileErr_fun01_sec06;
         else
            errSC = def_EOF_seqST;
      } /*If: gz file, might have more entries*/

      goto noErr_fun01_sec06;

   noErr_fun01_sec06:;
      seqSTPtr->seqLenSL = rmWhite_ulCp(seqSTPtr->seqStr);

      if(*typeSCPtr & def_fqType_gzSeqST)
      { /*If: fastq file*/
         seqSTPtr->qLenSL = rmWhite_ulCp(seqSTPtr->qStr);

         if(seqSTPtr->qLenSL < seqSTPtr->seqLenSL)
            goto badLine_fun01_sec06;
      } /*If: fastq file*/

      goto ret_fun01_sec06;

   memErr_fun01_sec06:;
      errSC = def_memErr_seqST;
      goto ret_fun01_sec06;

   fileErr_fun01_sec06:;
      errSC = def_fileErr_seqST;
      goto ret_fun01_sec06;

   badLine_fun01_sec06:;
      errSC = def_badLine_seqST;
      errSC |= def_fileErr_seqST;
      goto ret_fun01_sec06;

   ret_fun01_sec06:;
      return errSC;
} /*get_gzSeqST*/

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
:   FOR A PARTICSLAR PURPOSE AND NONINFRINGEMENT.  IN NO
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
:   FOR A PARTICSLAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
:   AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
:   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
