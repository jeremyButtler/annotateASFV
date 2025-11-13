/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' demux SOF: Start Of File
'   - has functions to demultiplexes reads by barcode
'   o header:
'     - guards and forward declarations
'   o fun01: sortBarcodes_demux
'     - sorts a the mapped barcodes by their barcode and
'       then start coordiante
'   o .c fun01: sortBarcodes_demux
'     - sorts the mapped barcodes by their barcode and
'       then by start coordiante
'   o .c fun02: sortStart_demux
'     - sorts the mapped barcodes by starting coordinates
'   o fun03: barcodeCoords_demux
'     - demux a read
'   o fun04: read_demux
'     - convert barcode coordinates to demuxed reads
'   o fun05: primer_demux
'     - get primer target regins from the input sequence
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - guards and forward declarations
\-------------------------------------------------------*/

#ifndef DE_MULTIPLEX_READS_H
#define DE_MULTIPLEX_READS_H

struct seqST;
struct tblST_kmerFind;
struct refST_kmerFind;

/*-------------------------------------------------------\
| Fun03: barcodeCoords_demux
|   - find the coordinates of barcodes in a read
| Input:
|   - lenSIPtr:
|     o signed int pointer to get length of returned array
|   - seqSTPtr:
|     o seqSTPtr struct pointer with read to demux
|   - barAryST:
|     o refST_kmerFind struct array with barcodes to demux
|   - barLenSI:
|     o number of barcodes in demux
|   - tblSTPtr:
|     o tblST_kmerFind struct pointer to use
|   - minPercScoreF:
|     o float with minimum percent score to keep a primer
|       mapping or discard
|     o recomend 0.9 (90%) to 0.95 (95%)
|     o do not go lower then 80% (otherwise to many
|       mappings kept)
|   - minScoreSL:
|     o signed long with minimum score to keep a primer
|       mapping or discard
|     o disable by using 0
|   - alnSetPtr:
|     o alnSet struct pointer with alignment settings
| Output:
|   - Modifies:
|     o lenSIPtr
|       - length of the returned array
|       - -1 if had two different barcodes overlapping
|       - -2 for memory errors
|   - Returns:
|     o array with the barcode index, start, end, and
|       score of each barcode
|       * index (n % 4 = 0) is barcode index
|       * index (n % 4 = 1) is the barcode start
|       * index (n % 4 = 2) is the barcode end
|       * index (n % 4 = 3) is the score
|       * ex: [index_1, start_1, end_1, score_1,
|              index_2, start_2 end_2, score_2,
|              index_3, start_3, end_3, score_3,
|              ...
|             ]
|     o for memory errors or no sequences (check lenSI)
\-------------------------------------------------------*/
signed int *
barcodeCoords_demux(
   signed int *lenSIPtr,    /*gets returned array length*/
   struct seqST *seqSTPtr,          /*has read to demux*/
   struct refST_kmerFind *barSTPtr, /*has barcodes*/
   signed int barLenSI,             /*number barcodes*/
   struct tblST_kmerFind *tblSTPtr, /*for searching*/
   float minPercScoreF,             /*min % score*/
   signed long minScoreSL,          /*min score to keep*/
   struct alnSet *alnSetPtr         /*settings*/
);

/*-------------------------------------------------------\
| Fun04: read_demux
|   - convert barcode coordinates to demuxed reads
| Input:
|   - seqSTPtr:
|     o seqSTPtr struct pointer with read to demux
|   - maxSplitSI:
|     o maximum number of splits before tossing read
|   - minDistSI:
|     o minimum distance between barcodes
|   - trimBl:
|     o 1: trim barcodes
|     o 0: do not trim barcodes
|   - coordsArySI:
|     o signed int array retured from barcodeCoords_demux
|       * index (n % 4 = 0) is barcode index
|       * index (n % 4 = 1) is the barcode start
|       * index (n % 4 = 2) is the barcode end
|       * index (n % 4 = 3) is the score
|   - coordLenSI:
|     o length of coordsArySI (number barcodes << 2)
|   - barFileStrAry:
|     o c-string array with output file names
|     o you need to check and make sure these files are
|       writeable before hand (I assume here)
| Output:
|   - Prints:
|     o demuxed and if requested split sequence to the
|       correct file in barFileStrAry
|   - Returns:
|     o 0 for no errors
|     o 1 if read had no barcodes
|     o 2 if read had to many barcodes
|     o 3 if read had barcodes at both ends; no easy split
|     o 4 if read had to barcodes right next to each other
\-------------------------------------------------------*/
signed int
read_demux(
   struct seqST *seqSTPtr,   /*read to demux*/
   signed int maxSplitSI,    /*max splits before discard*/
   signed int minDistSI, /*min distance between barcodes*/
   signed char trimBl,       /*1: trim barcodes*/
   signed int coordArySI[],  /*has barcode mappings*/
   signed int coordLenSI,    /*length of coordLenSI*/
   signed char *barFileStrAry[] /*output file names*/
);

/*-------------------------------------------------------\
| Fun05: primer_demux
|   - get primer target regins from the input sequence
| Input:
|   - seqSTPtr:
|     o seqSTPtr struct pointer with read to demux
|   - minDistSI:
|     o minimum distance between primers
|   - maxDistSI:
|     o maximum distance between primers
|   - coordsArySI:
|     o signed int array retured from barcodeCoords_demux
|       * index (n % 4 = 0) is barcode index
|       * index (n % 4 = 1) is the barcode start
|       * index (n % 4 = 2) is the barcode end
|       * index (n % 4 = 3) is the score
|         + a negative score means a reverse mapping
|   - coordLenSI:
|     o length of coordsArySI (number barcodes << 2)
|   - barSTPtr:
|     o refs_kmerFind struct pionter with the primer
|       ids and the index of their mates (if paired)
|   - outFILE:
|     o FILE pionter to print reads to
| Output:
|   - Prints:
|     o amplicons to outFILE
|   - Returns:
|     o number of amplicons found
|     o 0 if no amplicons
\-------------------------------------------------------*/
signed int
primer_demux(
   struct seqST *seqSTPtr,/*read to split into amplicons*/
   signed int minDistSI,  /*min distance between primers*/
   signed int maxDistSI,  /*max distance between primers*/
   signed int coordArySI[],/*has barcode mappings*/
   signed int coordLenSI, /*length of coordLenSI*/
   struct refST_kmerFind *barSTPtr,/*primer ids & mates*/
   void *outFILE          /*print sequences to*/
);

/*-------------------------------------------------------\
| Fun06: pGeneCoord_demux
|   - get coordinates of all genes found in target
| Input:
|   - seqSTPtr:
|     o seqSTPtr struct pointer with read to find genes
|       in
|   - coordsArySI:
|     o signed int array retured from barcodeCoords_demux
|       * index (n % 4 = 0) is barcode index
|       * index (n % 4 = 1) is the barcode start
|       * index (n % 4 = 2) is the barcode end
|       * index (n % 4 = 3) is the score
|         + a negative score means a reverse mapping
|   - coordLenSI:
|     o length of coordsArySI (number barcodes << 2)
|   - geneSTPtr:
|     o refST_kmerFind struct pionter with the genes
|   - headBlPtr:
|     o 1: print header and then set to 0
|     o 0: do not print the header
|   - geneSTPtr
|     o refST_kmerFind struct pionter with the genes
|       to find
|   - outFILE:
|     o FILE pionter to print gene coordinates
| Output:
|   - Prints:
|     o gene coordinates to outFILE
|   - Returns:
|     o number of amplicons found
|     o 0 if no genes
\-------------------------------------------------------*/
signed int
pGeneCoord_demux(
   struct seqST *seqSTPtr,/*read to split into amplicons*/
   signed int coordArySI[],/*has barcode mappings*/
   signed int coordLenSI, /*length of coordLenSI*/
   signed char *headBlPtr,/*1: print header + set to 0*/
   struct refST_kmerFind *geneSTPtr,/*genes searched*/
   void *outFILE          /*print sequences to*/
);

#endif

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconveint / not possible, this code is under the
:   MIT license
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
