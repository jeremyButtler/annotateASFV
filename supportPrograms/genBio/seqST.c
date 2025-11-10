/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' seqST SOF: Start Of File
'   - holds functions for reading in or manipulating
'     sequences
'   o header:
'     - included libraries
'   o .h st01: seqST
'     - holds an single sequence (fasta/fastq)
'   o fun02 getFq_seqST:
'     - reads a fastq sequence from a fastq file
'   o fun03 getFa_seqST:
'     - grabs the next read in the fasta file
'   o fun05 revComp_seqST:
'     - reverse complement a sequence
'   o fun06 blank_seqST:
'     - sets values in seqST to zero
'   O fun07 init_seqST:
'     - sets values in seqST to blank values
'   o fun08 freeStack_seqST:
'     - frees variables in an seqST (calls blank_seqST)
'   o fun09 freeHeap_seqST:
'     - frees an seqST strucuter (calls fredSeqSTStack)
'   o fun10: freeHeapAry_seqST
'     - frees an array of seqST's
'   o fun11 cpIdEndPad_seqST:
'      - copies read id to a buffer and adds in endIdC to
'        the end. If needed, this function will add right
'        padding of spaces to the end.
'   o fun12: cp_seqST
'     - copies an seqST structure
'   o fun13: swap_seqST
'     - swaps values in two seqST structs
'   o fun14: sort_seqST
'     - sorts an array of seqST structs
'   o fun15: search_seqST
'     - searchs for sequence in seqST struct array
'   o fun16: mkAry_seqST
'     - makes a seqST struct array
'   o fun17: realloc_seqST
'     - add more memory to a seqST struct array
'   o fun18: readFaFile_seqST
'     - get all sequences from a fasta file
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
#include "seqST.h"
#include "../genLib/ulCp.h"
#include "../genLib/fileFun.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - .h  #include "../genLib/endLine.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

/*-------------------------------------------------------\
| Fun02: getFq_seqST
|  - grabs next read in a fastq file
| Input:
|  - fqFILE:
|    o pointer to FILE handle to a fastq file to get the
|      next read from
|  - seqSTPtr:
|    o pointer to seqSTPtr structure to hold next read
| Output:
|  - Modifies:
|    o refStruct to hold the read in fastq entry & sets
|      its pointers
|  - Returns:
|     o 0: for no errors
|     o def_EOF_seqST: if EOF
|     o def_fileErr_seqST If file was not a fastq file
|     o def_badLine_seqST | def_fileErr_seqST for  a
|       invalid fastq entry
|     o def_memErr_seqST If malloc failed to find memory
\-------------------------------------------------------*/
signed char
getFq_seqST(
  void *fqFILE,           /*fastq file with sequence*/
  struct seqST *seqSTPtr  /*will hold one fastq entry*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun02 TOC: readRefFqSeq
   '  -  grabs the next read in the fastq file
   '  o fun02 sec01:
   '    - variable declarations
   '  o fun02 sec02:
   '    - read in header
   '  o fun02 sec03:
   '    - read in sequence
   '  o fun02 sec04:
   '    - read in quality score entry
   '  o fun02 sec05:
   '    - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun02 Sec01:
    ^  - variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    unsigned char errSC = 0;
    signed long lenSL = 0;
    signed long tmpSL = 0;
    signed long bytesSL = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun02 Sec02:
    ^  - read in the header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(! fqFILE)
        goto fileErr_fun02_sec05;

    if(! seqSTPtr->idStr)
    { /*If: need memory*/
       seqSTPtr->idStr =
          malloc(136 * sizeof(signed char));
       if(! seqSTPtr->idStr)
          goto memErr_fun02_sec05;
       seqSTPtr->idSizeSL = 128;
    } /*If: need memory*/

    seqSTPtr->idLenSL = 0;

    lenSL =
        getFullLine_fileFun(
            fqFILE,
            &seqSTPtr->idStr,
            &seqSTPtr->idSizeSL,
            &tmpSL,
            0
        ); /*get fastq header*/
    if(lenSL < 0)
       goto memErr_fun02_sec05;
    if(! lenSL)
       goto eof_fun02_sec05;

    while(lenSL && seqSTPtr->idStr[lenSL - 1] < 33)
       --lenSL;

    if(! lenSL)
       goto fileErr_fun02_sec05;/*blank line not allowed*/

    seqSTPtr->idLenSL = lenSL - 1;
    cpLen_ulCp(
       seqSTPtr->idStr,
       &seqSTPtr->idStr[1],
       seqSTPtr->idLenSL
    ); /*remove header symbol*/

    seqSTPtr->idStr[seqSTPtr->idLenSL] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 1] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 2] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 3] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 4] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 5] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 6] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 7] = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun02 Sec03:
    ^  - read in the sequence & spacer
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(! seqSTPtr->seqStr)
    { /*If: need memory*/
       seqSTPtr->seqStr =
          malloc(1032 * sizeof(signed char));
       if(! seqSTPtr->seqStr)
          goto memErr_fun02_sec05;
       seqSTPtr->seqSizeSL = 1024;

       seqSTPtr->seqStr[1024] = 0;
       seqSTPtr->seqStr[1025] = 0;
       seqSTPtr->seqStr[1026] = 0;
       seqSTPtr->seqStr[1027] = 0;
       seqSTPtr->seqStr[1028] = 0;
       seqSTPtr->seqStr[1029] = 0;
       seqSTPtr->seqStr[1030] = 0;
       seqSTPtr->seqStr[1031] = 0;
    } /*If: need memory*/


    seqSTPtr->seqLenSL = 0;
    bytesSL = 0;
    lenSL = 0;
    tmpSL = 0;

    do{ /*Loop: get sequence entry*/
       seqSTPtr->seqLenSL += lenSL;
       bytesSL += tmpSL;

       lenSL =
           getFullLine_fileFun(
               fqFILE,
               &seqSTPtr->seqStr,
               &seqSTPtr->seqSizeSL,
               &tmpSL,
               (signed long) seqSTPtr->seqLenSL
           ); /*get fastq header*/

       if(lenSL < 0)
          goto memErr_fun02_sec05;
       else if(! lenSL)
          goto lineErr_fun02_sec05;
    } while(seqSTPtr->seqStr[seqSTPtr->seqLenSL] != '+');
      /*Loop: get sequence entry*/

    seqSTPtr->seqStr[seqSTPtr->seqLenSL] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 1] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 2] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 3] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 4] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 5] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 6] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 7] = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun02 Sec04:
    ^  - read in q-score entry
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    seqSTPtr->qLenSL = 0;

    if(seqSTPtr->qSizeSL < seqSTPtr->seqLenSL)
    { /*If: need more memory*/
       if(seqSTPtr->qStr)
           free(seqSTPtr->qStr);
       seqSTPtr->qStr = 0;
       seqSTPtr->qSizeSL = 0;
    } /*If: need more memory*/

    if(! seqSTPtr->qStr)
    { /*If: assining memory*/
       seqSTPtr->qStr =
          malloc(
            (seqSTPtr->seqLenSL + 8) * sizeof(signed char)
          );
       if(! seqSTPtr->qStr)
          goto memErr_fun02_sec05;
       seqSTPtr->qSizeSL = seqSTPtr->seqLenSL;

       /*this is to aviod valgrind complates about
       `  unitialized values
       */
       seqSTPtr->qStr[seqSTPtr->seqLenSL] = 0;
       seqSTPtr->qStr[seqSTPtr->seqLenSL + 1] = 0;
       seqSTPtr->qStr[seqSTPtr->seqLenSL + 2] = 0;
       seqSTPtr->qStr[seqSTPtr->seqLenSL + 3] = 0;
       seqSTPtr->qStr[seqSTPtr->seqLenSL + 4] = 0;
       seqSTPtr->qStr[seqSTPtr->seqLenSL + 5] = 0;
       seqSTPtr->qStr[seqSTPtr->seqLenSL + 6] = 0;
       seqSTPtr->qStr[seqSTPtr->seqLenSL + 7] = 0;
    } /*If: assining memory*/

    seqSTPtr->qLenSL =
       fread(
          seqSTPtr->qStr,
          sizeof(signed char),
          bytesSL,
          (FILE *) fqFILE
       );


    seqSTPtr->qLenSL = rmWhite_ulCp(seqSTPtr->qStr);
    seqSTPtr->seqLenSL = rmWhite_ulCp(seqSTPtr->seqStr);
    seqSTPtr->offsetSL = 0;
    seqSTPtr->endAlnSL = seqSTPtr->seqLenSL;

    if(seqSTPtr->qLenSL < seqSTPtr->seqLenSL)
       goto lineErr_fun02_sec05;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun02 Sec05:
    ^   - return result
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    errSC = 0;
    goto ret_fun02_sec05;

    eof_fun02_sec05:;
       errSC = def_EOF_seqST;
       goto ret_fun02_sec05;

    memErr_fun02_sec05:;
       errSC = def_memErr_seqST;
       goto ret_fun02_sec05;

    fileErr_fun02_sec05:;
       errSC = def_fileErr_seqST;
       goto ret_fun02_sec05;

    lineErr_fun02_sec05:;
       errSC = def_fileErr_seqST;
       errSC |= def_badLine_seqST;
       goto ret_fun02_sec05;

    ret_fun02_sec05:;
       return errSC;
} /*getFq_seqST*/

/*-------------------------------------------------------\
| Fun03: getFa_seqST
|  -  grabs next read in fasta file
| Input:
|  - faFILE:
|    o pointer to FILE handle to a fasta file to get the
|      next read from
|  - seqSTPtr:
|    o pointer to seqSTPtr structure to hold next read
| Output:
|  - Modifies:
|    o seqSTPtr to hold one fasta entry
|  - Returns:
|     o 0: for no errors
|     o def_EOF_seqST: if EOF
|     o def_fileErr_seqST If file was not a fastq file
|     o def_badLine_seqST | def_fileErr_seqST for  a
|       invalid fasta entry
|     o def_memErr_seqST If malloc failed to find memory
| Note:
|   - This will remove new lines from the sequence.
|   - This will only remove spaces or tabs at the end of
|     each sequence entry, so "atg atc \n" will got to
|     "atcatc".
\-------------------------------------------------------*/
signed char
getFa_seqST(
  void *faFILE,           /*fasta file with sequence*/
  struct seqST *seqSTPtr  /*will hold one fastq entry*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun03 TOC: getFa_seqST
   '  - grabs the next read in the fasta file
   '  o fun03 sec01:
   '    - variable declarations
   '  o fun03 sec02:
   '    - read in header
   '  o fun03 sec03:
   '    - read in sequence
   '  o fun03 sec04:
   '    - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun03 Sec01:
    ^  - variable declarations
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    signed char tmpSC = 'C';
    unsigned char errSC = 0;
    signed long lenSL = 0;
    signed long tmpSL = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun03 Sec02:
    ^  - read in header
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(! faFILE)
       goto fileErr_fun03_sec04;

    if(! seqSTPtr->idStr)
    { /*If: need memory*/
       seqSTPtr->idStr =
          malloc(136 * sizeof(signed char));
       if(! seqSTPtr->idStr)
          goto memErr_fun03_sec04;
       seqSTPtr->idSizeSL = 128;
    } /*If: need memory*/

    seqSTPtr->idLenSL = 0;
    seqSTPtr->idStr[0] = 0;

    while(seqSTPtr->idStr[0] != '>')
    { /*Loop: find first header*/
       lenSL =
           getFullLine_fileFun(
               faFILE,
               &seqSTPtr->idStr,
               &seqSTPtr->idSizeSL,
               &tmpSL,
               0
           ); /*get fastq header*/

       if(! lenSL)
          goto eof_fun03_sec04;
       else if(lenSL < 0)
          goto memErr_fun03_sec04;

       tmpSL = 0;
       while(seqSTPtr->idStr[tmpSL] < 33 && tmpSL < lenSL)
          ++tmpSL;

       if(tmpSL == lenSL)
          continue; /*blank line*/
       else if(seqSTPtr->idStr[tmpSL] != '>')
          goto lineErr_fun03_sec04;
       else
          break;
    } /*Loop: find first header*/


    while(lenSL && seqSTPtr->idStr[lenSL - 1] < 33)
       --lenSL;
    seqSTPtr->idLenSL = lenSL;

    cpLen_ulCp(
       seqSTPtr->idStr,
       &seqSTPtr->idStr[1],
       seqSTPtr->idLenSL - 1
    ); /*remove header symbol*/

    seqSTPtr->idStr[seqSTPtr->idLenSL] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 1] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 2] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 3] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 4] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 5] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 6] = 0;
    seqSTPtr->idStr[seqSTPtr->idLenSL + 7] = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun03 Sec03:
    ^   - read in sequence
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    if(! seqSTPtr->seqStr)
    { /*If: need memory*/
       seqSTPtr->seqStr =
          malloc(1032 * sizeof(signed char));
       if(! seqSTPtr->seqStr)
          goto memErr_fun03_sec04;
       seqSTPtr->seqSizeSL = 1024;
    } /*If: need memory*/

    seqSTPtr->seqLenSL = 0;
    lenSL = 0;
    tmpSL = 0;

    do{ /*Loop: get sequence entry*/
       lenSL =
           getFullLine_fileFun(
               faFILE,
               &seqSTPtr->seqStr,
               &seqSTPtr->seqSizeSL,
               &tmpSL,
               (signed long) seqSTPtr->seqLenSL
           ); /*get fastq header*/

       if(lenSL < 0)
          goto memErr_fun03_sec04;
       else if(! lenSL && seqSTPtr->seqLenSL)
          goto done_fun03_sec04;
       else if(! lenSL)
          goto lineErr_fun03_sec04;

       seqSTPtr->seqLenSL += lenSL;
       tmpSC = fgetc(faFILE); /*get next valid char*/
       while( ((unsigned char) tmpSC) < 33 )
          tmpSC = fgetc(faFILE); /*get next valid char*/

       if(tmpSC == '>')
          ungetc(tmpSC, faFILE);
       else if(tmpSC == EOF)
          goto done_fun03_sec04;
       else
          ungetc(tmpSC, faFILE);
    } while(tmpSC != '>');
      /*Loop: get sequence entry*/

    seqSTPtr->seqStr[seqSTPtr->seqLenSL] = 0;

    /*this is so always have eight initialized values for
    `  ulCp
    */
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 1] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 2] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 3] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 4] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 5] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 6] = 0;
    seqSTPtr->seqStr[seqSTPtr->seqLenSL + 7] = 0;

    /*<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<\
    ^ Fun03 Sec04:
    ^   - check for errors and return
    \>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>*/

    done_fun03_sec04:;
       seqSTPtr->seqLenSL =
          rmWhite_ulCp(seqSTPtr->seqStr);

       seqSTPtr->offsetSL = 0;
       seqSTPtr->endAlnSL = seqSTPtr->seqLenSL - 1;

       errSC = 0;
       goto ret_fun03_sec04;

    eof_fun03_sec04:;
       if(! seqSTPtr->seqStr) ;
       else if(! *seqSTPtr->seqStr) ;
       else
          seqSTPtr->seqLenSL =
             rmWhite_ulCp(seqSTPtr->seqStr);

       seqSTPtr->offsetSL = 0;
       seqSTPtr->endAlnSL = seqSTPtr->seqLenSL;

       errSC = def_EOF_seqST;
       goto ret_fun03_sec04;

    memErr_fun03_sec04:;
       errSC = def_memErr_seqST;
       goto ret_fun03_sec04;

    fileErr_fun03_sec04:;
       errSC = def_fileErr_seqST;
       goto ret_fun03_sec04;

    lineErr_fun03_sec04:;
       errSC = def_fileErr_seqST;
       errSC = def_badLine_seqST;
       goto ret_fun03_sec04;

    ret_fun03_sec04:;
       return errSC;
} /*getFa_seqST*/

/*-------------------------------------------------------\
| Fun05: revComp_seqST
|  - reverse complement a sequence
| Input:
|  - seqSTPtr:
|    o pointer to seqSTPtr structure with a sequence to
|      reverse complement
| Output:
|  - Modfies
|    o reverse complements the sequence entry and if
|      reverses the Q-score entry
\-------------------------------------------------------*/
void
revComp_seqST(
  struct seqST *seqSTPtr  /*sequence to reverse comp*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
    ' SOF: Start Of File
    '   o fun05 sec01:
    '     - variable declarations
    '   o fun05 sec02:
    '     - check if have q-score entry
    '   o fun05 sec03:
    '     - reverse complement sequence
    '   o fun05 sec04:
    '     - reverse complement final base and return
    \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun05 Sec01:
    ^   - variable declarations
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    signed char *seqStr = seqSTPtr->seqStr;

    signed char *endStr =
       seqSTPtr->seqStr + seqSTPtr->seqLenSL - 1;

    signed char *qStr = 0;
    signed char *qEndStr = 0;

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun05 Sec02:
    ^   - check if have q-score entry
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    if(seqSTPtr->qStr != 0 && *seqSTPtr->qStr != '\0')
    { /*If have a Q-score entry*/
       qStr = seqSTPtr->qStr;
       qEndStr = seqSTPtr->qStr + seqSTPtr->qLenSL - 1;
    } /*If have a Q-score entry*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun05 Sec03:
    ^   - reverse complement sequence
    ^   o fun05 sec03 sub01:
    ^     - reverse complement start base + start loop
    ^   o fun05 sec03 sub02:
    ^     - reverse complement end base
    ^   o fun05 sec03 sub03:
    ^     - swap start and end bases
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /****************************************************\
    * Fun05 Sec03 Sub01:
    *   - reverse complement start base + start loop
    \****************************************************/

    while(endStr > seqStr)
    { /*While I have bases to reverse complement*/
       switch(*seqStr & (~32)) /*Make sure upper case*/
       { /*switch, reverse complement*/
           case 'A': *seqStr = 'T'; break;
           case 'C': *seqStr = 'G'; break;
           case 'G': *seqStr = 'C'; break;
           case 'T': *seqStr = 'A'; break;
           case 'U': *seqStr = 'A'; break;
           case 'W': *seqStr = 'W'; break;
           case 'S': *seqStr = 'S'; break;
           case 'M': *seqStr = 'K'; break;  /*A/C->T/G*/
           case 'K': *seqStr = 'M'; break;  /*T/G->A/C*/
           case 'R': *seqStr = 'Y'; break;  /*A/G->T/C*/
           case 'Y': *seqStr = 'R'; break;  /*T/C->A/G*/
           case 'B': *seqStr = 'V'; break;/*C/G/T->G/C/A*/
           case 'D': *seqStr = 'H'; break;/*G/T/A->C/A/T*/
           case 'H': *seqStr = 'D'; break;/*C/A/T->G/T/A*/
           case 'V': *seqStr = 'B'; break;/*A/C/G->T/G/C*/
           case 'N': *seqStr = 'N'; break;/*A/C/G/T*/
           default: *seqStr = 0;
       } /*switch, reverse complement*/

       /*************************************************\
       * Fun05 Sec03 Sub02:
       *   - reverse complement end base
       \*************************************************/

       switch(*endStr & (~32)) /*Make sure upper case*/
       { /*switch, reverse complement*/
           case 'A': *endStr = 'T'; break;
           case 'C': *endStr = 'G'; break;
           case 'G': *endStr = 'C'; break;
           case 'T': *endStr = 'A'; break;
           case 'U': *endStr = 'A'; break;
           case 'W': *endStr = 'W'; break;
           case 'S': *endStr = 'S'; break;
           case 'M': *endStr = 'K'; break;  /*A/C->T/G*/
           case 'K': *endStr = 'M'; break;  /*T/G->A/C*/
           case 'R': *endStr = 'Y'; break;  /*A/G->T/C*/
           case 'Y': *endStr = 'R'; break;  /*T/C->A/G*/
           case 'B': *endStr = 'V'; break;/*C/G/T->G/C/A*/
           case 'D': *endStr = 'H'; break;/*G/T/A->C/A/T*/
           case 'H': *endStr = 'D'; break;/*C/A/T->G/T/A*/
           case 'V': *endStr = 'B'; break;/*A/C/G->T/G/C*/
           case 'N': *endStr = 'N'; break;/*A/C/G/T*/
           default: *endStr = 0;
       } /*switch, reverse complement*/

       /*************************************************\
       * Fun05 Sec03 Sub03:
       *   - swap start and end bases
       \*************************************************/

       *seqStr ^= *endStr;
       *endStr ^= *seqStr;
       *seqStr ^= *endStr;
       
       if(qStr)
       { /*If I also need to swap Q-scores*/
         *qStr ^= *qEndStr;
         *qEndStr ^= *qStr;
         *qStr ^= *qEndStr;

         ++qStr;
         --qEndStr;
       } /*If I also need to swap Q-scores*/
       
       ++seqStr;
       --endStr;
    } /*While I have bases to reverse complement*/

    /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
    ^ Fun05 Sec04:
    ^   - reverse complement final base and return
    \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

    /*Check if ended on same base; if so complement base*/
    if(endStr == seqStr)
    { /*If: ending on same base (not swaped)*/
       switch(*seqStr & (~32)) /*Make sure upper case*/
       { /*switch, reverse complement*/
           case 'A': *seqStr = 'T'; break;
           case 'C': *seqStr = 'G'; break;
           case 'G': *seqStr = 'C'; break;
           case 'T': *seqStr = 'A'; break;
           case 'U': *seqStr = 'A'; break;
           case 'W': *seqStr = 'W'; break;
           case 'S': *seqStr = 'S'; break;
           case 'M': *seqStr = 'K'; break;  /*A/C->T/G*/
           case 'K': *seqStr = 'M'; break;  /*T/G->A/C*/
           case 'R': *seqStr = 'Y'; break;  /*A/G->T/C*/
           case 'Y': *seqStr = 'R'; break;  /*T/C->A/G*/
           case 'B': *seqStr = 'V'; break;/*C/G/T->G/C/A*/
           case 'D': *seqStr = 'H'; break;/*G/T/A->C/A/T*/
           case 'H': *seqStr = 'D'; break;/*C/A/T->G/T/A*/
           case 'V': *seqStr = 'B'; break;/*A/C/G->T/G/C*/
           case 'N': *seqStr = 'N'; break;/*A/C/G/T*/
           default: *seqStr = 0;
       } /*switch, reverse complement*/
    } /*If: ending on same base (not swaped)*/

    return;
} /*revComp_seqST*/

/*-------------------------------------------------------\
| Fun06: blank_seqST
|  - sets values in seqST to blank values
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer to blank
| Output:
|  - Modifies
|    o sets id, sequence, and q-score entries to start
|      with '\0' and the id, sequence, and q-score lengths
|      to 0. This does not change the buffer lengths.
\-------------------------------------------------------*/
void
blank_seqST(
   struct seqST *seqSTPtr
){
   if(seqSTPtr->idStr)
      *seqSTPtr->idStr = '\0';

   seqSTPtr->idLenSL = 0;

   if(seqSTPtr->seqStr)
      *seqSTPtr->seqStr='\0';

   seqSTPtr->seqLenSL = 0;

   if(seqSTPtr->qStr)
      *seqSTPtr->qStr = '\0';

   seqSTPtr->qLenSL = 0;

   seqSTPtr->offsetSL = 0;
   seqSTPtr->endAlnSL = 0;
} /*blank_seqST*/

/*-------------------------------------------------------\
| Fun07: init_seqST
|  - sets vlues in seqST to zero
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer to initialize
| Output:
|  - Modifies
|    o all values in seqST to be 0
\-------------------------------------------------------*/
void
init_seqST(
   struct seqST *seqSTPtr
){
   seqSTPtr->idStr = 0;
   seqSTPtr->idSizeSL = 0;

   seqSTPtr->seqStr = 0;
   seqSTPtr->seqSizeSL = 0;

   seqSTPtr->qStr = 0;
   seqSTPtr->qSizeSL = 0;

   blank_seqST(seqSTPtr);
} /*init_seqST*/

/*-------------------------------------------------------\
| Fun08: freeStack_seqST
|  - frees the variables in an seqST strucuter
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer that has variables
|      to free
| Output:
|  - Frees
|    o seqSTPtr-idStr and sets to 0
|    o seqSTPtr-seqStr and sets to 0
|    o seqSTPtr-qStr and sets to 0
\-------------------------------------------------------*/
void
freeStack_seqST(
   struct seqST *seqSTPtr
){
   if(seqSTPtr)
   { /*If: I have something to free*/
      free(seqSTPtr->idStr);
      seqSTPtr->idStr = 0;

      free(seqSTPtr->seqStr);
      seqSTPtr->seqStr = 0;

      free(seqSTPtr->qStr);
      seqSTPtr->qStr = 0;

      blank_seqST(seqSTPtr);
   } /*If: I have something to free*/
} /*freeStack_seqST*/

/*-------------------------------------------------------\
| Fun09: freeHeap_seqST
|  - frees the seqST strucuter
| Input:
|  - seqSTPtr:
|    o pointer to an seqST structuer to free
| Output:
|  - frees
|    o seqST from memory
\-------------------------------------------------------*/
void
freeHeap_seqST(
   struct seqST *seqSTPtr
){
   freeStack_seqST(seqSTPtr); /*checks if have null*/
   free(seqSTPtr);
} /*freeHeap_seqST*/

/*-------------------------------------------------------\
| Fun10: freeHeapAry_seqST
|  - frees an array of seqST's
| Input:
|  - seqAryST:
|    o pointer to an arrya of seqST structuer to free
|  - numSeqSI:
|    o number of seqSTs in the seqST array
|     (seqSTAryPtr)
| Output:
|  - Frees
|    o seqST from memory
\-------------------------------------------------------*/
void
freeHeapAry_seqST(
   struct seqST *seqAryST,
   signed long numSeqSL
){
   signed long slFree = 0;

   for(
      slFree = 0;
      slFree < numSeqSL;
      ++slFree
   ){ /*Loop: Free the individual structures*/
      freeStack_seqST(&seqAryST[slFree]);
   } /*Loop: Free the individual structures*/

   free(seqAryST);
} /*freeHeapAry_seqST*/

/*-------------------------------------------------------\
| Fun11: cpIdEndPad_seqST
|  - copies read id to a buffer and adds in endIdC to
|    the end. If needed, this function will add right
|    padding of spaces to the end.
| Input:
|  - seqSTPtr:
|    o pointer to seqST with read id to copy
|  - buffStr:
|    o c-string to copy read id to
|  - endIdC:
|    o character to add in to mark the end of the read id
|  - padRI:
|    o amount of padding to add to the right of the read
|      id
| Output:
|  - Modifies:
|    o buffStr to hold the sequence id + endIdC. If the
|      id + endIdC is shorter than padRI, the copied id is
|      padded with spaces on the right till it is the same
|      size as padRI.
|  - Returns:
|    o pointer to end of copied id or padding if padding
|      was applied.
\-------------------------------------------------------*/
signed char *
cpIdEndPad_seqST(
   struct seqST *seqSTPtr,/*has read id to copy*/
   signed char *buffStr, /*buffer to add read id to*/
   signed char endIdC,   /*padding; id end (0 to skip)*/
   signed int padRI      /*padding; id start*/
){
   signed int idCntI = 0;
   signed char *tmpStr = seqSTPtr->idStr;

   if(*tmpStr == '>') ++tmpStr; /*Get of header*/

   while(*tmpStr != '\0')
   { /*While I have a read id to copy over*/
      *buffStr = *tmpStr;
      ++buffStr;
      ++tmpStr;
      ++idCntI;
   } /*While I have a read id to copy over*/

   if(endIdC)
   { /*If adding in a character after the read id*/
      *buffStr = endIdC;
      ++buffStr; /*If not a null*/
      ++idCntI;
   } /*If adding in a character after the read id*/

   while(idCntI < padRI)
   { /*While I have padding to add to the end*/
      *buffStr = ' ';
      ++buffStr;
      ++idCntI;
   } /*While I have padding to add to the end*/

   return buffStr;
} /*cpIdEndPad_seqST*/

/*-------------------------------------------------------\
| Fun12: cp_seqST
|   - copies an seqST structure
| Input:
|   - dupSeqST:
|     o pointer to an seqST to copy to (duplicate)
|   - cpSeqST:
|     o pointer to an seqST to copy
| Output:
|   - Modifies:
|     o dupSeqST to be an copy of cp_seqST
|   - returns:
|     o 0 for no errors
|     o def_memErr_seqST for an memory error
\-------------------------------------------------------*/
signed char
cp_seqST(
   struct seqST *dupSeqST, /*copy to (duplicate)*/
   struct seqST *cpSeqST   /*seqST to copy*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun12 TOC:
   '   - copies an seqST structure
   '   o fun12 sec01:
   '     - blank and copy aligment settings
   '   o fun12 sec02:
   '     - copy q-score entry
   '   o fun12 sec03:
   '     - copy sequence entry
   '   o fun12 sec04:
   '     - copy read id
   '   o fun12 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec01:
   ^   - blank and copy aligment settings
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   blank_seqST(dupSeqST);
   dupSeqST->offsetSL = cpSeqST->offsetSL;
   dupSeqST->endAlnSL = cpSeqST->endAlnSL;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec02:
   ^   - copy q-score entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/
   
   if(cpSeqST->qStr)
   { /*If: have q-score entry*/
      if(
         ! dupSeqST->qStr
         || dupSeqST->qSizeSL < cpSeqST->qLenSL
      ){ /*If: have to resize buffer*/

         if(dupSeqST->qStr)
         { /*If: need to resize buffer*/
            free(dupSeqST->qStr);
            dupSeqST->qStr = 0;
         } /*If: need to resize buffer*/

         dupSeqST->qSizeSL = cpSeqST->qLenSL;

         dupSeqST->qStr =
            malloc(
                 (dupSeqST->qSizeSL + 1)
               * sizeof(signed char)
            );

         if(! dupSeqST)
            goto memErr_fun12;
      } /*If: have to resize buffer*/
      
      dupSeqST->qLenSL = cpSeqST->qLenSL;

      cpLen_ulCp(
         dupSeqST->qStr,
         cpSeqST->qStr,
         dupSeqST->seqLenSL
      );
   } /*If: have q-score entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec03:
   ^   - copy sequence entry
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
      ! dupSeqST->seqStr
      || dupSeqST->seqSizeSL < cpSeqST->seqLenSL
   ){ /*If: have to resize buffer*/

      if(dupSeqST->seqStr)
      { /*If: need to resize buffer*/
         free(dupSeqST->seqStr);
         dupSeqST->seqStr = 0;
      } /*If: need to resize buffer*/

      dupSeqST->seqSizeSL = cpSeqST->seqLenSL;

      dupSeqST->seqStr =
         malloc(
              (dupSeqST->seqSizeSL + 1)
            * sizeof(signed char)
         );

      if(! dupSeqST)
         goto memErr_fun12;
   } /*If: have to resize buffer*/

   dupSeqST->seqLenSL = cpSeqST->seqLenSL;

   cpLen_ulCp(
      dupSeqST->seqStr,
      cpSeqST->seqStr,
      dupSeqST->seqLenSL
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec04:
   ^   - copy read id
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(
      ! dupSeqST->idStr
      || dupSeqST->idSizeSL < cpSeqST->idLenSL
   ){ /*If: have to resize buffer*/

      if(dupSeqST->idStr)
      { /*If: need to resize buffer*/
         free(dupSeqST->idStr);
         dupSeqST->idStr = 0;
      } /*If: need to resize buffer*/

      dupSeqST->idSizeSL = cpSeqST->idLenSL;

      dupSeqST->idStr =
         malloc(
              (dupSeqST->idSizeSL + 1)
            * sizeof(signed char)
         );

      if(! dupSeqST)
         goto memErr_fun12;
   } /*If: have to resize buffer*/
   
   dupSeqST->idLenSL = cpSeqST->idLenSL;

   cpLen_ulCp(
      dupSeqST->idStr,
      cpSeqST->idStr,
      dupSeqST->idLenSL
   );

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun12 Sec05:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   return 0;

   memErr_fun12:;
   return def_memErr_seqST;
} /*cp_seqST*/

/*-------------------------------------------------------\
| Fun13: swap_seqST
|   - swaps values in two seqST structs
| Input:
|   - firstSTPtr:
|     o pointer first seqST struct to swap
|   - secSTPtr:
|     o pointer second seqST struct to swap
| Output:
|   - Modifies:
|     o firstSTPtr to have values from secSTPtr
|     o secSTPtr to have values from firstSTPtr
\-------------------------------------------------------*/
void
swap_seqST(
   struct seqST *firstSTPtr,
   struct seqST *secSTPtr
){
   signed char *swapStr = 0;

   /*read id*/
   swapStr = firstSTPtr->idStr;
   firstSTPtr->idStr = secSTPtr->idStr;
   secSTPtr->idStr = swapStr;

   firstSTPtr->idLenSL ^= secSTPtr->idLenSL;
   secSTPtr->idLenSL ^= firstSTPtr->idLenSL;
   firstSTPtr->idLenSL ^= secSTPtr->idLenSL;

   firstSTPtr->idSizeSL ^= secSTPtr->idSizeSL;
   secSTPtr->idSizeSL ^= firstSTPtr->idSizeSL;
   firstSTPtr->idSizeSL ^= secSTPtr->idSizeSL;


   /*sequence*/
   swapStr = firstSTPtr->seqStr;
   firstSTPtr->seqStr = secSTPtr->seqStr;
   secSTPtr->seqStr = swapStr;

   firstSTPtr->seqLenSL ^= secSTPtr->seqLenSL;
   secSTPtr->seqLenSL ^= firstSTPtr->seqLenSL;
   firstSTPtr->seqLenSL ^= secSTPtr->seqLenSL;

   firstSTPtr->seqSizeSL ^= secSTPtr->seqSizeSL;
   secSTPtr->seqSizeSL ^= firstSTPtr->seqSizeSL;
   firstSTPtr->seqSizeSL ^= secSTPtr->seqSizeSL;


   /*q-score entry*/
   swapStr = firstSTPtr->qStr;
   firstSTPtr->qStr = secSTPtr->qStr;
   secSTPtr->qStr = swapStr;

   firstSTPtr->qLenSL ^= secSTPtr->qLenSL;
   secSTPtr->qLenSL ^= firstSTPtr->qLenSL;
   firstSTPtr->qLenSL ^= secSTPtr->qLenSL;

   firstSTPtr->qSizeSL ^= secSTPtr->qSizeSL;
   secSTPtr->qSizeSL ^= firstSTPtr->qSizeSL;
   firstSTPtr->qSizeSL ^= secSTPtr->qSizeSL;


   /*alignent variables*/
   firstSTPtr->offsetSL ^= secSTPtr->offsetSL;
   secSTPtr->offsetSL ^= firstSTPtr->offsetSL;
   firstSTPtr->offsetSL ^= secSTPtr->offsetSL;

   firstSTPtr->endAlnSL ^= secSTPtr->endAlnSL;
   secSTPtr->endAlnSL ^= firstSTPtr->endAlnSL;
   firstSTPtr->endAlnSL ^= secSTPtr->endAlnSL;
} /*swap_seqST*/

/*-------------------------------------------------------\
| Fun14: sort_seqST
|   - sorts an array of seqST structs
| Input:
|   - arySTPtr:
|     o pointer to array of seqST struct to sort
|   - lenSL:
|     o length to sort in array
| Output:
|   - Modifies:
|     o arySTPtr to be sorted by read/ref id
\-------------------------------------------------------*/
void
sort_seqST(
   struct seqST *arySTPtr, /*array to sort*/
   signed long lenSL       /*length of array*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun14 TOC:
   '   - sorts an array of seqST structs
   '   o fun14 sec01:
   '     - variable declerations
   '   o fun14 sec02:
   '     - find the number of rounds to sort for
   '   o fun14 sec03:
   '     - sort the arrays
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec01:
   ^   - variable declerations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed long cmpSL = 0;

   /*Number of sorting rounds*/
   signed long subSL = 0;
   signed long nextSL = 0;
   signed long lastSL = 0;
   signed long onSL = 0;

   /*Variables to incurment loops*/
   signed long slIndex = 0;
   signed long slElm = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec02:
   ^   - find the max search value (number rounds to sort)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(! lenSL)
      return;

   /*Recursion formsla: h[0] = 1, h[n] = 3 * h[n - 1] +1*/
   subSL = 1; /*Initialzie first array*/

   while(subSL < lenSL - 1)
      subSL = (3 * subSL) + 1;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun14 Sec03:
   ^   - sort arrays
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(subSL > 0)
   { /*Loop: all rounds*/
      for(
         slIndex = 0;
         slIndex <= subSL;
         ++slIndex
      ){ /*Loop: though sub array*/
         slElm = slIndex;

         for(
            slElm = slIndex;
            slElm + subSL < lenSL;
            slElm += subSL
         ){ /*Loop: swap elements in subarray*/
            nextSL = slElm + subSL;

            cmpSL =
               eqlNull_ulCp(
                  arySTPtr[slElm].idStr,
                  arySTPtr[nextSL].idStr
               ); /*See if need to sort structs*/

            if(cmpSL > 0)
            { /*If I need to swap an element*/
               swap_seqST(
                  &arySTPtr[slElm],
                  &arySTPtr[nextSL]
               );

               lastSL = slElm;
               onSL = slElm;

               while(lastSL >= subSL)
               { /*loop: move swapped element back*/
                  lastSL -= subSL;

                  cmpSL =
                     eqlNull_ulCp(
                        arySTPtr[onSL].idStr,
                        arySTPtr[lastSL].idStr
                     ); /*See if need to sort structs*/

                  if(cmpSL > 0)
                     break; /*Positioned the element*/
   
                  swap_seqST(
                     &arySTPtr[onSL],
                     &arySTPtr[lastSL]
                  );

                  onSL = lastSL;
               } /*Loop: move swapped element back*/
            } /*If I need to swap elements*/
         } /*Loop: swap elements in subarray*/
      } /*Loop: though sub array*/

      subSL = (subSL - 1) / 3; /*Move to next round*/
   } /*Loop: all rounds*/
} /*sort_seqST*/

/*-------------------------------------------------------\
| Fun15: search_seqST
|  - searchs for sequence in seqST struct array
| Input:
|  - seqArySTPtr:
|    o pionter to seqST stuct array to search
|  - qryStr:
|    o query to find
|  - lenSL:
|    o length of seqArySTPtr (index 1)
| Output:
|  - Returns:
|    o index of qryStr in seqArySTPtr
|    o -1 if qryUL is not in ulAry
\-------------------------------------------------------*/
signed long
search_seqST(
   struct seqST *seqArySTPtr, /*sequence array to search*/
   signed char *qryStr,       /*query to hunt for*/
   signed long lenSL          /*length to search in*/
){
   signed long cmpSL = 0;
   signed long midSL = 0;
   signed long rightSL = lenSL - 1;
   signed long leftSL = 0;

   while(leftSL <= rightSL)
   { /*Loop: Search for the querys index*/
      midSL = (leftSL + rightSL) >> 1;

      cmpSL =
         eqlNull_ulCp(
            qryStr,
            seqArySTPtr[midSL].idStr
         ); /*compare query to sequence at midpiont*/

      if(cmpSL > 0)
         leftSL = midSL + 1;

      else if(cmpSL < 0)
         rightSL = midSL - 1;

      else
        return midSL;
   } /*Loop: Search for the querys index*/

   return -1; /*query not found*/
} /*search_seqST*/

/*-------------------------------------------------------\
| Fun16: mkAry_seqST
|  - makes a seqST struct array
| Input:
|  - sizeSL:
|    o size of seqST struct array
| Output:
|  - Returns:
|    o pointer to seqST struct array
|    o 0 for memory errors
\-------------------------------------------------------*/
struct seqST *
mkAry_seqST(
   signed long sizeSL         /*new length (size)*/
){
   struct seqST *retHeapAryST = 0;

   retHeapAryST = malloc(sizeSL * sizeof(struct seqST));

   if(! retHeapAryST)
      goto memErr_fun16;

   while(sizeSL > 0)
   { /*Loop: initialize structs*/
      --sizeSL;   /*start at index 1*/
      init_seqST( &retHeapAryST[sizeSL] );
   } /*Loop: initialize structs*/

   return retHeapAryST;

   memErr_fun16:;
      return 0;
} /*mkAry_seqST*/

/*-------------------------------------------------------\
| Fun17: realloc_seqST
|  - add more memory to a seqST struct array
| Input:
|  - seqArySTPtr:
|    o pionter to seqST stuct array to expand
|  - lenSL:
|    o number of items in array (index 1)
|    o this marks the start of initialization
|  - sizeSL:
|    o new size of seqST struct array
| Output:
|  - Returns:
|    o 0 for no errors
|    o def_memErr_seqST for memory errors
\-------------------------------------------------------*/
signed char
realloc_seqST(
   struct seqST **seqArySTPtr,/*sequence array to expand*/
   signed long lenSL,         /*number items in array*/
   signed long sizeSL         /*new length (size)*/
){
   struct seqST *tmpSTPtr = 0;

   tmpSTPtr =
      realloc(
         *seqArySTPtr,
         sizeSL * sizeof(struct seqST)
      );

   if(! tmpSTPtr)
      goto memErr_fun16;

   *seqArySTPtr = tmpSTPtr;

   while(lenSL < sizeSL)
   { /*Loop: initialize unused arrays*/
      init_seqST( &(*seqArySTPtr)[lenSL] );
      ++lenSL;
   } /*Loop: initialize unused arrays*/

   return 0;

   memErr_fun16:;
      return def_memErr_seqST;
} /*realloc_seqST*/

/*-------------------------------------------------------\
| Fun18: readFaFile_seqST
|   - get all sequences from a fasta file
| Input:
|   - faFILE:
|     o FILE pointer to file to get sequences from
|   - lenSLPtr:
|     o pointer to signed long to hold number sequences
|       in returned array (also in faFILE)
|   - sizeSLPtr:
|     o pointer to signed long to hold returned array size
|   - errSCPtr:
|     o pointer to signed char to hold errors
| Output:
|   - Modifies:
|     o faFILE to be at end of file
|     o lenSLPtr to have number sequences in faFILE/array
|     o sizeSLPtr to have size of returned array
|     o errSCPtr to have errors
|       * 0 for no errors
|       * def_fileErr_seqST if had file error
|       * def_badLine_seqST | def_fileErr_seqST for
|         invalid fasta entry
|       * def_memErr_seqST for memory errors
|   - Returns:
|     o seqST struct array with sequences
|     o 0 for errors
\-------------------------------------------------------*/
struct seqST *
readFaFile_seqST(
   void *faFILE,           /*file to get reads from*/
   signed long *lenSLPtr,  /*number items in array*/
   signed long *sizeSLPtr, /*length of array*/
   signed char *errSCPtr   /*holds error values*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun18 TOC:
   '   - get all sequences from a fasta file
   '   o fun18 sec01:
   '     - variable declarations
   '   o fun18 sec02:
   '     - allocate memory and get first sequence
   '   o fun18 sec03:
   '     - read rest of fasta file
   '   o fun18 sec04:
   '     - return results
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   struct seqST *retHeapAryST = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec02:
   ^   - allocate memory and get first sequence
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retHeapAryST = mkAry_seqST(16);

   if(! retHeapAryST)
      goto memErr_fun18_sec04;

   *sizeSLPtr = 16;

   *errSCPtr =
      getFa_seqST((FILE *) faFILE, &retHeapAryST[0]);

   if(*errSCPtr == def_EOF_seqST)
      goto fileErr_fun18_sec04; /*no sequences in file*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec03:
   ^   - read rest of fasta file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *lenSLPtr = 1;
   
   while(! *errSCPtr)
   { /*Loop: read in references*/

      if(*lenSLPtr >= *sizeSLPtr)
      { /*If: need more memory*/
         *errSCPtr =
            realloc_seqST(
               &retHeapAryST,
               *lenSLPtr,
               (*sizeSLPtr << 1)
            );

         if(*errSCPtr)
            goto memErr_fun18_sec04;

         *sizeSLPtr <<= 1;
      } /*If: need more memory*/

      *errSCPtr =
         getFa_seqST(
            (FILE *) faFILE,
            &retHeapAryST[*lenSLPtr]
         );

      if(*errSCPtr == def_EOF_seqST)
      { /*If: hit end of file*/
         freeStack_seqST(&retHeapAryST[*lenSLPtr]);
            /*make sure no loose ends*/
         break; /*finished*/
      } /*If: hit end of file*/

      ++(*lenSLPtr);
   } /*Loop: read in references*/


   if(*errSCPtr != def_EOF_seqST)
      goto errClean_fun18_sec04;
      /*error reading file; already corect value*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec04:
   ^   - return results
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *errSCPtr = 0;
   goto ret_fun18_sec04;

   memErr_fun18_sec04:;
      *errSCPtr = def_memErr_seqST;
      goto errClean_fun18_sec04;

   fileErr_fun18_sec04:;
      *errSCPtr = def_fileErr_seqST;
      goto errClean_fun18_sec04;

   errClean_fun18_sec04:;
      if(retHeapAryST)
         freeHeapAry_seqST(
            retHeapAryST,
            *sizeSLPtr
         );
      retHeapAryST = 0;
      goto ret_fun18_sec04;

   ret_fun18_sec04:;
      return retHeapAryST;
} /*readFaFile_seqST*/

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
