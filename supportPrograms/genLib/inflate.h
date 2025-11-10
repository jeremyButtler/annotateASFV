/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' inflate SOF: Start Of File
'   - inflate a fastx.gz file (not best, but should work?)
'   o header:
'     - included libraries and tables
'   o .h st01: tree_unizp
'     - has tree to decode ziped data
'   o fun01: blank_tree_inflate
'     - blank values in a tree_inflate struct
'   o fun02: init_tree_inflate
'     - initializes values in a tree_inflate struct
'   o .c fun03: freeStack_tree_inflate
'     - would free stack, but tree_inflate has no payload
'     - not needed, but habit
'   o .c fun04: freeHeap_tree_inflate
'     - frees tree_inflate struct; does not free children
'     - not needed, but habit
'   o fun05: freeTreeHeap_tree_inflate
'     - frees an entire tree_inflate struct tree
'   o fun06: insert_tree_inflate
'     - inserts new huffman node into tree_inflate struct
'   o fun07: blank_file_inflate
'     - blanks a file_inflate structure
'   o fun08: init_file_inflate
'     - initializes a file_inflate structure
'   o fun09: freeStack_file_inflate
'     - frees variabls in a file_inflate structure
'   o fun10: freeHeap_file_inflate
'     - frees a file_inflate structure
'   o fun13: getZipBuff_file_inflate
'     - checks and then if needed, reads new data in
'       compressed data buffer
'   o .c fun14: getBit_inflate
'     - get next bit in array and incument to next
'       position
'   o .c fun15: getNBits_inflate
'     - get next bits in array and incurment position
'   o .c fun16: getBytes_inflate
'     - get next bytes in array & incument to next pos
'   o .c fun17: getValue_inflate
'     - gets a value from zlib stream and converts to a
'       value using a huffman tree
'   o fun18: dictToTree_inflate
'     - converts a dictionary to a huffman tree
'   o fun19: getTree_inflate
'     - gets huffman trees for inflation
'   o fun20: getFixedTree_inflate
'     - gets huffman fixed trees for fixed inflation
'   o fun21: block_inflate
'     - inflates data in a block
'   o fun22: cpBytes_inflate
'     - copies bytes from file_inflate struct a buffer
'   o fun23: get_inflate
'     - inflate blocks until have required input
'   o fun24: getc_file_inflate
'     - get one charater or one line break, which can two
'       be up to two characters from file
'   o fun25: unget_file_inflate
'     - add one character back to the buffer; only call
'       once per read
'   o fun26: next_file_inflate
'     - move to next deflate/gz entry in the file
'   o fun27: add_file_inflate
'     - decompresses zip data
'   o license:
'     - licensing for this code (public domain / mit)
' Credit:
'   - I realied on a blog pretty heavily to get this code
'     https://pyokagan.name/blog/2019-10-18-zlibinflate/
'   - good C guide for zlib and gzip
`     https://www.infinitepartitions.com/art001.html
'   - might be good for compression
'     https://www.euccas.me/zlib/
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - defined variables and header guards
|   o header sec01:
|     - error values and header gaurd
|   o header sec02:
|     - data/file type variables and maximum copy sizes
\-------------------------------------------------------*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec01:
^   - error values and header gaurd
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#ifndef INFLATE_GZ_H
#define INFLATE_GZ_H

#define def_memErr_inflate 1
#define def_badSymbol_inflate 2
#define def_eofErr_inflate 3   /*premature EOF*/
#define def_eof_inflate 4      /*valid EOF*/
#define def_endBlock_inflate 5
#define def_fullBuff_inflate 6 /*buffer is full*/
#define def_badBlock_inflate 7 
#define def_fileErr_inflate 8
#define def_badType_inflate 9
#define def_badWindow_inflate 10
#define def_badCheck_inflate 11
#define def_dictErr_inflate 12

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Header Sec02:
^   - data/file type variables and maximum copy sizes
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define def_noBlock_inflate 4 /*emtpy value*/
   /*values 0 to 2 are reserved for zlib*/

/*file types*/
#define def_zlib_inflate 0 /*not checked*/
#define def_gz_inflate 1

#define def_maxLen_inflate 258
   /*maximum copy length for deflate*/
#define def_maxDist_inflate 32768
   /*32768 is the maximum call back length for inflate*/
#define def_window_inflate (def_maxDist_inflate + 2)
   /*my window size, the + 2 is to allow `\n` or `\r\n`*/

/*-------------------------------------------------------\
| ST01: tree_inflate
|   - has tree to decode ziped data
\-------------------------------------------------------*/
typedef struct tree_inflate
{
   unsigned short valUS;
   struct tree_inflate *zeroST;
   struct tree_inflate *oneST;
}tree_inflate;

/*-------------------------------------------------------\
| ST02: file_inflate
|   - has inflate buffers and trees
\-------------------------------------------------------*/
typedef struct file_inflate
{
   signed char finalBl;     /*1: is final block*/
   signed char blockFlagSC; /*holds block type on*/
   unsigned short bytesLeftUS; /*bytes left in block*/
   signed char eofBl;       /*marks if at end of file*/

   /*has compressed bytes to uncompress; zip buffer*/
   unsigned char zipByteUC; /*byte to process*/
   unsigned char zipMaskUC; /*has bit to extract*/

   /*has uncompressed bytes to copy to users buffer*/
   signed char buffStr[def_window_inflate + 8];
   unsigned short buffSizeUS;/*number bytes in buffer*/
   unsigned short buffLenUS;/*uncompress position*/
   unsigned short buffPosUS;/*last copied byte*/

   unsigned char checkSumAryUC[4];/*for finding checksum*/

   struct tree_inflate *lenTreeST; /*length huffman tree*/
   struct tree_inflate *distTreeST;/*distance huffman*/

   void *zipFILE; /*file to get data from*/


   /*these are for gz format*/
   unsigned char timeAryUC[4];   /*time stamp*/
   signed char nameStr[256];     /*file name*/
   unsigned char nameLenUC;      /*length of file name*/
   signed char typeSC;           /*file type*/
      /*0 is deflate; 1 is gzip*/
   signed char binBl;            /*1 is binary; 0 text*/

   unsigned int crc32UI; /*for crc32 checksums*/
   unsigned int lenUI;   /*bytes uncompressed*/
}file_inflate;

/*-------------------------------------------------------\
| Fun01: blank_tree_inflate
|   - blank values in a tree_inflate struct
| Input:
|   - treeSTPtr:
|      o pointer to tree_inflate struct to blank
| Output:
|   - Modifies:
|     o valSS in treeSTPtr to be 0
\-------------------------------------------------------*/
void
blank_tree_inflate(
   struct tree_inflate *treeSTPtr
);

/*-------------------------------------------------------\
| Fun02: init_tree_inflate
|   - initializes values in a tree_inflate struct
| Input:
|   - treeSTPtr:
|      o pointer to tree_inflate struct to initialize
| Output:
|   - Modifies:
|     o valSS in treeSTPtr to be 0
|     o leftST and rightST in treeSTPtr to be 0
\-------------------------------------------------------*/
void
init_tree_inflate(
   struct tree_inflate *treeSTPtr
);

/*-------------------------------------------------------\
| Fun05: freeTreeHeap_tree_inflate
|   - frees an entire tree_inflate struct tree
| Input:
|   - treeSTPtr:
|      o pointer to tree_inflate struct free
| Output:
|   - Frees:
|     o treeSTPtr and all its children (you must set
|       to null)
\-------------------------------------------------------*/
void
freeTreeHeap_tree_inflate(
   struct tree_inflate *treeSTPtr
);

/*-------------------------------------------------------\
| Fun06: insert_tree_inflate
|   - inserts new huffman node into a tree_inflate struct
| Input:
|   - treeSTPtr:
|      o pointer to tree_inflate struct insert into
|   - charSC:
|     o character chode represents
|   - codeSS:
|     o code for symbol adding
|   - lenCodeSS:
|     o number bits in codeSS
| Output:
|   - Modifies:
|     o treeSTPtr to have new symbol in it (somewere)
|   - Returns:
|     o 0 for no errors
|     o def_memErr_inflate for memory errors
\-------------------------------------------------------*/
signed char
insert_tree_inflate(
   struct tree_inflate *treeSTPtr, /*tree to insert in*/
   signed short valSS,           /*value to add*/
   signed short codeSS,          /*code for valSS*/
   signed short lenCodeSS        /*bits in code*/
);

/*-------------------------------------------------------\
| Fun07: blank_file_inflate
|   - blanks a file_inflate structure
| Input:
|   - fileSTPtr:
|     o file_inflate structure pointer to blank
|   - closeFileBl:
|     o 2: close zip file and set to 0
|     o 1: seek start of zip file
|       * stdin input is always closed (set to 0)
|     o 0: do nothing
| Output:
|   - Modifies:
|      o zipLenSL, zipPosSL, in fileSTPtr to 0
|      o buffLenUS, and buffPosUS in fileSTPtr to 0
|      o finalBl in fileSTPtr to 0
|      o blockFlagSC in fileSTPtr to def_noBlock_inflate
|      o zipBuffStr in fileSTPtr to be '\0' (if exist)
|      o zipFILE in fileSTPtr to be closed and set to 0
\-------------------------------------------------------*/
void
blank_file_inflate(
   struct file_inflate *fileSTPtr,
   signed char closeFileBl
);

/*-------------------------------------------------------\
| Fun08: init_file_inflate
|   - initializes a file_inflate structure
| Input:
|   - fileSTPtr:
|     o file_inflate structure pointer to initialize
| Output:
|   - Modifies:
|      o all values in fileSTPtr to be 0/null
\-------------------------------------------------------*/
void
init_file_inflate(
   struct file_inflate *fileSTPtr
);

/*-------------------------------------------------------\
| Fun09: freeStack_file_inflate
|   - frees variabls in a file_inflate structure
| Input:
|   - fileSTPtr:
|     o file_inflate structure pointer with vars to free
| Output:
|   - Modifies:
|      o frees and sets all variables and files to 0
\-------------------------------------------------------*/
void
freeStack_file_inflate(
   struct file_inflate *fileSTPtr
);

/*-------------------------------------------------------\
| Fun10: freeHeap_file_inflate
|   - frees a file_inflate structure
| Input:
|   - fileSTPtr:
|     o file_inflate structure pointer to free
| Output:
|   - Modifies:
|      o frees fileSTPtr, but you must set to 0/null
|        * includes files in fileSTPtr
\-------------------------------------------------------*/
void
freeHeap_file_inflate(
   struct file_inflate *fileSTPtr
);

/*-------------------------------------------------------\
| Fun13: getZipBuff_file_inflate
|   - checks and then if needed, reads new data in
|     compressed data buffer
| Input:
|   - fileSTPtr:
|     o file_inflate struct pointer to add compressed
|       data to
|   - offsetSL:
|     o number of bytes needed for next read
|     o forces read before reached size limit; max is 1024
| Output:
|   - Modifies:
|     o zipAryUC in fileSTPtr to have new data
|     o zipLenSL in fileSTPtr to have new data length
|     o zipPosSL in fileSTPtr to have new position
|     o zipFILE in fileSTPtr to be on next set of bytes
|   - Returns:
|     o 0 for no erorrs
|     o def_eof_inflate for end of file (EOF)
\-------------------------------------------------------*/
signed char
getZipBuff_file_inflate(
   struct file_inflate *fileSTPtr,
   signed long offsetSL /*old data to keep*/
);

/*-------------------------------------------------------\
| Fun18: dictToTree_inflate
|   - converts a dictionary to a huffman tree
| Input:
|   - lenArySS:
|     o unsigned short array of bit lengths of each code
|   - lenSS:
|     o number codes/code lengths (dictionary length)
| Output:
|   - Returns:
|     o tree_inflate pointer to huffman tree
|     o 0 for memory errors
\-------------------------------------------------------*/
struct tree_inflate *
dictToTree_inflate(
   unsigned short *lenAryUS,       /*length of each code*/
   unsigned short lenUS            /*number of codes*/
);

/*-------------------------------------------------------\
| Fun19: getTree_inflate
|   - gets huffman trees for inflation
| Input:
|   - fileSTPtr:
|     o file_inflate struct pointer with zip buffers/file
|       and to add tree to
| Output:
|   - Modifies:
|     o zipAryUC, zipPosSL, zipLenSL an to be on the
|       next positions
|     o distTreeST and lenTreeST to have huffman trees
|   - Returns:
|     o 0 for no errors
|     o def_memErr_inflate for memory errors
|     o def_badSymbol_inflate for invalid tree symbol
|     o def_eofErr_inflate for early EOF
\-------------------------------------------------------*/
signed char
getTree_inflate(
   struct file_inflate *fileSTPtr
);

/*-------------------------------------------------------\
| Fun20: getFixedTree_inflate
|   - gets huffman fixed trees for fixed inflation
| Input:
|   - fileSTPtr:
|     o file_inflate struct with huffman trees
| Output:
|   - Modifies:
|     o distTreeST and lenTreeST in fileSTPtr to have the
|       distance and length trees
|   - Returns:
|     o 0 for no errors
|     o def_memErr_inflate for memory errors
\-------------------------------------------------------*/
signed char
getFixedTree_inflate(
   struct file_inflate *fileSTPtr
);

/*-------------------------------------------------------\
| Fun21: block_inflate
|   - inflates data in a block
| Input:
|   - fileSTPtr:
|     o file_inflate struct pointer with zip buffer/file
|       to read
|   - errSCPtr:
|     o signed char pointer to hold errors
| Output:
|   - Modifies:
|     o zipPosSL, zipLenSL, zipAryUC, and zipFILE in
|       fileSTPtr to be set to/have next set of compressed
|       bytes
|       * reads in new data and resets if needed
|     o buffStr and buffLenUS in fileSTPtr to have
|       uncompressed data in block
|     o buffPosUS in fileSTPtr to be on corrected position
|       if had to make more room in buffer
|     o errSCPtr to hold error type
|       * 0 for no errors
|       * def_endBlock_inflate if at end of block
|       * def_eofErr_inflate for early EOF
|       * def_fullBuff_inflate if buffer is full; means
|         you need to copy bytes from buffer and incurment
|         buffPosUS
|   - Returns:
|     o number of converted bytes
\-------------------------------------------------------*/
signed long
block_inflate(
   struct file_inflate *fileSTPtr,/*has zip to inflate*/
   signed char *errSCPtr          /*gets errors*/
);

/*-------------------------------------------------------\
| Fun22: cpBytes_inflate
|   - copies bytes from file_inflate struct a buffer
| Input:
|   - fileSTPtr:
|     o file_inflate struct pointer with buffer to copy
|       bytes from
|   - outStr:
|     o c-string to copy bytes to
|   - posSL:
|     o offset (position) to start copy at in outStr
|   - newLineBl:
|     o 4: there are two extra bytes in buffer
|       * these are not included in bytesUS and are used
|         only for adding a line break to the end of a
|         buffer
|     o 1: copy till newline found
|     o 0: copy input number of bytes
|   - bytesUS:
|     o maximum number bytes to copy
| Output:
|   - Returns:
|     o number of bytes copied
|     o 0 if could no bytes in buffer to copy
| Note:
|   - this copies the line ending for the users system,
|     not the line ending stored in the zip file
\-------------------------------------------------------*/
unsigned short
cpBytes_inflate(
   struct file_inflate *fileSTPtr,
   signed char *outStr,    /*buffer to hold output*/
   signed long posSL,      /*offset at in outStr*/
   signed char newLineBl,  /*convert line breaks to OS*/
   unsigned short bytesUS  /*number bytes to get*/
);

/*-------------------------------------------------------\
| Fun23: get_inflate
|   - inflate blocks until have required input
| Input:
|   - fileSTPtr:
|     o file_inflate struct pointer to uncompres data with
|   - bytesSL:
|     o maximum bytes to get
|   - endBl:
|     o 8: skip buffer grab step (forces refill if no
|          buffPosSL is one off buffLenSL)
|     o 4: have two extra bytes at end of buffer for line
|          breaks (only used if line breaks)
|     o 2: copy until have full block
|     o 1: copy until new line: '\r', '\n', '\r\n', '\n\r'
|       or until bytesSL is reached
|     o 0: do not end at new lines, only bytesSL
|   - outStr:
|     o c-stirng to get uncompressed data
|   - errSCPtr
|     o signe char pointer to hold errors
| Output:
|   - Modifies:
|     o everything in fileSTPtr to be incurmented to end
|       of read in buffer
|     o outStr to have copied buffer contents
|     o errSCPtr
|       * 0 for no errors
|       * def_eof_inflate for end of final block
|         i. this is not file EOF, just decompress EOF
|       * def_memErr_inflate for memory errors
|       * def_badSymbol_inflate for bad symbol in zip file
|       * def_badBlock_inflate for unrecognized blocks
|       * def_eofErr_inflate for early EOF
|   - Returns:
|     o number bytes read into blockArySC
|     o < 0 for error (excluding def_eof_inflate)
\-------------------------------------------------------*/
signed long
get_inflate(
   struct file_inflate *fileSTPtr,
   signed long bytesSL,
   signed char endBl,
   signed char *outStr,
   signed char *errSCPtr
);

/*-------------------------------------------------------\
| Fun24: getc_file_inflate
|   - get one charater or one line break (two characters)
|     from file
| Input:
|   - fileSTPtr:
|     o file_inflate struct pointer with buffer to getc
|   - getArySC:
|     o c-string (size 3) to hold character or line break
|       with null at end
|       * for non-line breaks 1 byte + null byte = 2 bytes
|       * line break 1-2 bytes (OS) + null byte = 3 bytes
| Output:
|   - Modifies:
|     o buffStr in fileSTPtr to be moved one character or
|       one break
|     o buffPosSL in fileSTPtr to be one character
|       or break forward
|   - Returns:
|     o number characters in getArySC
|     o def_eof_inflate * -1 for end of final block
|       * this is not file EOF, just decompress EOF
|     o def_memErr_inflate * -1 for memory errors
|     o def_badSymbol_inflate * -1 for bad symbol in zip
|       file
|     o def_badBlock_inflate * -1 for unrecognized blocks
|     o def_eofErr_inflate * -1 for early EOF
\----------:--------------------------------------------*/
signed char
getc_file_inflate(
   struct file_inflate *fileSTPtr,
   signed char getArySC[]
);

/*-------------------------------------------------------\
| Fun24: getc_file_inflate
|   - get one charater or one line break (two characters)
|     from file
| Input:
|   - fileSTPtr:
|     o file_inflate struct pointer with buffer to getc
|   - getArySC:
|     o c-string (size 3) to hold character or line break
|       with null at end
|       * for non-line breaks 1 byte + null byte = 2 bytes
|       * line break 1-2 bytes (OS) + null byte = 3 bytes
| Output:
|   - Modifies:
|     o buffStr in fileSTPtr to be moved one character or
|       one break
|     o buffPosSL in fileSTPtr to be one character
|       or break forward
|   - Returns:
|     o number characters in getArySC
|     o def_eof_inflate * -1 for end of final block
|       * this is not file EOF, just decompress EOF
|     o def_memErr_inflate * -1 for memory errors
|     o def_badSymbol_inflate * -1 for bad symbol in zip
|       file
|     o def_badBlock_inflate * -1 for unrecognized blocks
|     o def_eofErr_inflate * -1 for early EOF
\----------:--------------------------------------------*/
signed char
getc_file_inflate(
   struct file_inflate *fileSTPtr,
   signed char getArySC[]
);

/*-------------------------------------------------------\
| Fun25: unget_file_inflate
|   - add one character back to the buffer; only call once
|     per read
| Input:
|   - fileSTPtr:
|     o file_inflate struct pointer with buffer to unget
|   - ungetSC:
|     o character to unget
| Output:
|   - Modifies:
|     o buffStr in fileSTPtr to have ungetSC
|     o buffPosSL in fileSTPtr to be one character back
\-------------------------------------------------------*/
void
unget_file_inflate(
   struct file_inflate *fileSTPtr,
   signed char ungetSC
);

/*-------------------------------------------------------\
| Fun26: next_file_inflate
|   - move to next deflate/gz entry in the file
| Input:
|   - fileSTPtr:
|     o file_inflate struct with gz file to get next file
|       for
|    - firstBytesArySC:
|      o signed char array of size two with the first two
|        bytes in the file
|        * this allows you to check for the gz header
|          before calling this function
|      o or 0/null for no bytes
|    - newBl:
|      o 1 or not 0: is new file
|      o 0: at end of a gz entry, move to next
| Output:
|   - Modifies:
|     o fileSTPtr to be on the next entry
|   - Returns:
|     o 0 for no errors
|     o def_eof_inflate if at ned of file
|     o def_badType_inflate for wrong compression type
|     o def_badWindow_inflate if window size to large
|     o def_badCheck_inflate if failed a checksum
|     o def_dictErr_inflate if zip uses dictionary
|     o def_eofErr_inflate if hit end of file to early
\-------------------------------------------------------*/
signed char
next_file_inflate(
   struct file_inflate *fileSTPtr, /*has next entry*/
   signed char *firstBytesArySC,   /*1st 2 bytes in file*/
   signed char newBl               /*1: first time*/
);

/*-------------------------------------------------------\
| Fun27: add_file_inflate
|   - adds a file to an file_inflate struct
| Input:
|   - zipFILE:
|     o FILE pointer to add in (set to 0/null after)
|       * always set this to 0/null after finshed, since
|         freeing the file_inflate struct will use fclose
|         (or for stdin, set to 0)
|       * 0/null will make file error
|   - fileSTPtr:
|     o file_inflate struct pointer to hold new file
|   - firstBytesArySC:
|     o signed char array of size two with the first two
|       bytes in the file
|       * this allows you to check for the gz header
|         before calling this function
|     o or 0/null for no bytes
| Output:
|   - Modifies:
|     o fileSTPtr to have zipFILE and header bytes of
|       zipFILE (zipPosSL, zipLenSL, zipAryUC)
|     o zipFILE to be on entries after header
|     o timeAryUC, nameStr, nameLenUC, typeSC, and binBl
|       to have their assigned data (if gz file supplied)
|   - Returns:
|     o 0 for no errors
|     o def_fileErr_inflate if no file input
|     o def_badType_inflate for wrong compression type
|     o def_badWindow_inflate if window size to large
|     o def_badCheck_inflate if failed header checksum
|     o def_dictErr_inflate if zip uses dictionary
|     o def_eofErr_inflate if hit end of file early, but
|       after first two bytes (something in file)
|     o def_eof_inflate if hit end of file on first two
|       bytes (likely blank file)
\-------------------------------------------------------*/
signed char
add_file_inflate(
   void *zipFILE,                  /*file to add*/
   struct file_inflate *fileSTPtr, /*gets file*/
   signed char *firstBytesArySC    /*1st 2 bytes in file*/
);

#endif

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
