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
'   - I realied on two blogs and one code base pretty
'     heavily to get this code
'     - https://pyokagan.name/blog/2019-10-18-zlibinflate/
`     - https://www.infinitepartitions.com/art001.html
'     - https://github.com/davecom/nflate
'       - only one with working huffman tree building step
'     - github copilot hepled with crc32 sums
'   - good C guide for zlib and gzip
'   - might be good for compression
'     - https://www.euccas.me/zlib/
'   - has huffman table, I would like to come back and
'     implement this to increase my speed
'     - https://commandlinefanatic.com/cgi-bin/showarticle.cgi?article=art007
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries and tables
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include "inflate.h"
#include <stdio.h>

#include "genMath.h"
#include "ulCp.h"
#include "endin.h"
#include "checkSum.h"

/*only use .h files*/
#include "endLine.h"

signed short lenExtraTbl_inflate[] =
   {
       11,  13,  15,  17,  19,  23, 27,
       31,  35,  43,  51,  59,  67, 83,
       99, 115, 131, 163, 195, 227 
   };

signed short distExtraTbl_inflate[] =
   {
           4,     6,     8,    12,    16,    24,
           32,   48,    64,    96,   128,   192,
          256,   384,  512,   768,  1024,  1536, 
         2048,  3072, 4096,  6144,  8192, 12288,
        16384, 24576 
   };

signed short codeLenTbl_inflate[] =
   {
      16, 17, 18,  0,  8,  7,  9,  6, 10, 5,
      11,  4, 12,  3, 13,  2, 14,  1, 15
   };

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
){
   if(! treeSTPtr)
     return;

   treeSTPtr->valUS = (unsigned short) -1;
} /*blank_tree_inflate*/

/*-------------------------------------------------------\
| Fun02: init_tree_inflate
|   - initializes values in a tree_inflate struct
| Input:
|   - treeSTPtr:
|      o pointer to tree_inflate struct to initialize
| Output:
|   - Modifies:
|     o valSS in treeSTPtr to be 0
|     o zeroST and oneST in treeSTPtr to be 0
\-------------------------------------------------------*/
void
init_tree_inflate(
   struct tree_inflate *treeSTPtr
){
   if(! treeSTPtr)
     return;

   treeSTPtr->oneST = 0;
   treeSTPtr->zeroST = 0;
   blank_tree_inflate(treeSTPtr);
} /*init_tree_inflate*/

/*-------------------------------------------------------\
| Fun03: freeStack_tree_inflate
|   - would free stack, but tree_inflate has no payload
| Input:
|   - treeSTPtr:
|      o pointer to tree_inflate struct stack free
| Output:
|   - Frees:
|     o initializes treeSTPtr (make sure children freeded)
\-------------------------------------------------------*/
void
freeStack_tree_inflate(
   struct tree_inflate *treeSTPtr
){
   if(! treeSTPtr)
      return;
   init_tree_inflate(treeSTPtr);
} /*freeStack_tree_inflate*/

/*-------------------------------------------------------\
| Fun04: freeHeap_tree_inflate
|   - frees a tree_inflate struct (does not free children)
| Input:
|   - treeSTPtr:
|      o pointer to tree_inflate struct free
| Output:
|   - Frees:
|     o treeSTPtr and no children (you must set to null)
\-------------------------------------------------------*/
void
freeHeap_tree_inflate(
   struct tree_inflate *treeSTPtr
){
   if(! treeSTPtr)
      return;
   freeStack_tree_inflate(treeSTPtr);
   free(treeSTPtr);
}

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun05 TOC:
   '   - frees an entire tree_inflate struct tree
   '   o fun05 sec01:
   '     - variable declarations and check if have input
   '   o fun05 sec02:
   '     - set up for loop (last tree and move down)
   '   o fun05 sec03:
   '     - free nodes
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec01:
   ^   - variable declarations and check if have input
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   struct tree_inflate *lastSTPtr = 0;
   struct tree_inflate *swapSTPtr = 0;

   if(! treeSTPtr)
     return;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec02:
   ^   - set up for loop (last tree and move down)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   lastSTPtr = treeSTPtr;

   if(treeSTPtr->zeroST)
      treeSTPtr = treeSTPtr->zeroST;
   else if(treeSTPtr->oneST)
      treeSTPtr = treeSTPtr->oneST;
   else
   { /*Else: only root node*/
      free(treeSTPtr);
      return;
   } /*Else: only root node*/

   lastSTPtr->zeroST = 0; /*so gets 0 at end*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun05 Sec03:
   ^   - free nodes
   ^ +-->* always set next nodes left child to parent node
   ^ |     o this allows me to work backwards
   ^ |   * always move left when possible
   ^ |   * when hit leaf (no childrend) free leaf, then
   ^ |     move to parent node using zeroST (left child)
   ^ |+->* if parent node has no right child, free
   ^ ||
   ^ |+--* repeat until parent has right child or no nodes
   ^ |   * if have node, move to its right child
   ^ +---* if have node; repeat until no nodes left
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   while(treeSTPtr)
   { /*Loop: free tree*/
      if(! treeSTPtr->zeroST)
      { /*If: could be at leaf*/
  
         if(! treeSTPtr->oneST)
         { /*If: at leaf*/
            free_fun05_sec03:;

            while(! treeSTPtr->oneST)
            { /*Loop: free node*/
               if(! lastSTPtr)
                  ;
               else if(lastSTPtr->oneST == treeSTPtr)
                  lastSTPtr->oneST = 0;

               freeHeap_tree_inflate(treeSTPtr);
               treeSTPtr = 0;
               treeSTPtr = lastSTPtr;

               if(lastSTPtr)
                  lastSTPtr = lastSTPtr->zeroST;
               else if(! treeSTPtr)
                  break;
            } /*Loop: free node*/

            if(! treeSTPtr)
               break;
            else
               continue;
         } /*If: at leaf*/

         else
         { /*Else: not at leaf*/
            treeSTPtr->zeroST = lastSTPtr;
            lastSTPtr = treeSTPtr;
            treeSTPtr = treeSTPtr->oneST;
         } /*Else: not at leaf*/
      } /*If: could be at leaf*/

      else if(! (treeSTPtr->zeroST == lastSTPtr) )
      { /*Else: need to swap*/
         swapSTPtr = treeSTPtr->zeroST;
         treeSTPtr->zeroST = lastSTPtr;
         lastSTPtr = treeSTPtr;
         treeSTPtr = swapSTPtr;
      } /*Else: need to swap*/

      else if(treeSTPtr->oneST)
      { /*Else If: parent already linked, move to oneST*/
         lastSTPtr = treeSTPtr;
         treeSTPtr = treeSTPtr->oneST;
      } /*Else If: parent already linked, move to oneST*/

      else
          goto free_fun05_sec03; /*at leaf*/
   } /*Loop: free tree*/
} /*freeTreeHeap_tree_inflate*/

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
|      o zipAryUC in fileSTPtr to be '\0' (if exist)
|      o zipFILE in fileSTPtr to be closed and set to 0
\-------------------------------------------------------*/
void
blank_file_inflate(
   struct file_inflate *fileSTPtr,
   signed char closeFileBl
){
   if(! fileSTPtr)
      return;

   fileSTPtr->finalBl = 0;
   fileSTPtr->blockFlagSC = def_noBlock_inflate;
   fileSTPtr->bytesLeftUS = 0;
   fileSTPtr->eofBl = 0;

   fileSTPtr->zipByteUC = 0;
   fileSTPtr->zipMaskUC = 8;

   fileSTPtr->buffStr[def_window_inflate] = '\0';
   fileSTPtr->buffStr[0] = '\0';
   fileSTPtr->buffSizeUS = 0;
   fileSTPtr->buffLenUS = 0;
   fileSTPtr->buffPosUS = 0;

   fileSTPtr->checkSumAryUC[0] = 0;
   fileSTPtr->checkSumAryUC[1] = 0;
   fileSTPtr->checkSumAryUC[2] = 0;
   fileSTPtr->checkSumAryUC[3] = 0;

   if(closeFileBl & 2)
   { /*If: closing the file*/
      if(
            fileSTPtr->zipFILE
         && fileSTPtr->zipFILE != stdout
         && fileSTPtr->zipFILE != stderr
         && fileSTPtr->zipFILE != stdin
      ) fclose((FILE *) fileSTPtr->zipFILE);
      fileSTPtr->zipFILE = 0;
   } /*If: closing the file*/

   else if(closeFileBl & 1)
   { /*Else If: seeking start of file*/
      if(
               fileSTPtr->zipFILE
            && fileSTPtr->zipFILE != stdout
            && fileSTPtr->zipFILE != stderr
            && fileSTPtr->zipFILE != stdin
      ) fseek(fileSTPtr->zipFILE, 0, SEEK_SET);

      else
         fileSTPtr->zipFILE = 0;
         /*stdin, can not rewind*/
   } /*Else If: seeking start of file*/


   fileSTPtr->nameStr[0] = '\0';
   fileSTPtr->nameLenUC = 0;

   fileSTPtr->typeSC = def_zlib_inflate;

   fileSTPtr->timeAryUC[0] = 0;
   fileSTPtr->timeAryUC[1] = 0;
   fileSTPtr->timeAryUC[2] = 0;
   fileSTPtr->timeAryUC[3] = 0;

   fileSTPtr->binBl = 0;
   fileSTPtr->crc32UI = (unsigned int) -1;
   fileSTPtr->lenUI = 0;
} /*blank_file_inflate*/

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
){
   if(! fileSTPtr)
      return;

   fileSTPtr->lenTreeST = 0;
   fileSTPtr->distTreeST = 0;
   fileSTPtr->zipFILE = 0;

   blank_file_inflate(fileSTPtr, 2);
} /*init_file_inflate*/

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
){
   if(! fileSTPtr)
      return;

   if(fileSTPtr->lenTreeST)
      freeTreeHeap_tree_inflate(fileSTPtr->lenTreeST);
   if(fileSTPtr->distTreeST)
      freeTreeHeap_tree_inflate(fileSTPtr->distTreeST);

   if(
         fileSTPtr->zipFILE
      && fileSTPtr->zipFILE != stdout
      && fileSTPtr->zipFILE != stderr
      && fileSTPtr->zipFILE != stdin
   ) fclose((FILE *) fileSTPtr->zipFILE);
   fileSTPtr->zipFILE = 0;

   init_file_inflate(fileSTPtr);
} /*freeStack_file_inflate*/

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
){
   if(! fileSTPtr)
      return;
   freeStack_file_inflate(fileSTPtr);
   free(fileSTPtr);
} /*freeHeap_file_inflate*/

/*-------------------------------------------------------\
| Fun14: getBit_inflate
|   - get next bit in array and incument to next position
| Input
|   - fileSTPtr:
|     o file_inflate struct pointer to get bit from
| Output:
|   - Modifies:
|     o eofBl in fileSTPtr to be 1 if at EOF
|       - also does nothing if eofBl is not 0
|     o sets zipMaskUC in fileSTPtr to next bit to extract
|     o sets zipFILE to next byte (if need new byte)
|   - Returns:
|     o next bit in aryMac or 0 if hit end of file
\-------------------------------------------------------*/
unsigned char
getBit_inflate(
   struct file_inflate *fileSTPtr
){
   unsigned char bitUC = 0;

   if(fileSTPtr->eofBl)
      return 0;

   if(fileSTPtr->zipMaskUC & 8)
   { /*If: need next byte*/
      fileSTPtr->zipMaskUC = 0;

      if(
         ! fread(
            &fileSTPtr->zipByteUC,
            sizeof(unsigned char),
            1,
            fileSTPtr->zipFILE
         )
      ) fileSTPtr->eofBl = 1;

      else
         fileSTPtr->eofBl = 0;
   } /*If: need next byte*/

   bitUC =
      (fileSTPtr->zipByteUC >> fileSTPtr->zipMaskUC) & 1;
   ++fileSTPtr->zipMaskUC;

   return bitUC;
} /*getBit_inflate*/

/*-------------------------------------------------------\
| Fun15: getNBits_inflate
|   - get next bits in array and incurment position
| Input
|   - fileSTPtr:
|     o file_inflate sturct pointer with buffer to get
|       bits from (and file to get more data from)
|   - bitsSI:
|     o number of bits to read
|   - errSCPtr:
|     o signed char pointer to hold any errors
| Output:
|   - Modifies:
|     o zipMaskUC in fileSTPtr to be set for the next bit
|     o zipByteUC to be next byte (if needed)
|     o zipFILE in fileSTPtr to be on next unused byte
|     o errSCPtr
|       * 0 for no errors
|       * def_eofErr_inflate for early EOF
|   - Returns:
|     o next byte in aryMac
| Note:
|   - I think zip will never ask for more than 16 bits
|     at a time (unless it is asking for bytes)
\-------------------------------------------------------*/
signed short
getNBits_inflate(
   struct file_inflate *fileSTPtr,
   signed int bitsSI,
   signed char *errSCPtr
){
   signed short retSS = 0;
   signed int siBit = 0;

   /*needs to be in little endin order*/
   for(siBit = 0; siBit < bitsSI; ++siBit)
   { /*Loop: get bits*/
      retSS |= ( getBit_inflate(fileSTPtr) << siBit );

      if(fileSTPtr->eofBl)
         goto eofErr_fun15;
   }  /*Loop: get bits*/

   *errSCPtr = 0;
   return retSS;

   eofErr_fun15:;
      *errSCPtr = def_eofErr_inflate;
      return retSS;
} /*getNBits_inflate*/

/*-------------------------------------------------------\
| Fun16: getBytes_inflate
|   - get next bytes in array & incument to next position
| Input
|   - fileSTPtr:
|     o file_inflate sturct pointer with buffer to get
|       bytes from (and file to get more data from)
|   - bytesSI:
|     o number of bytes to read
|   - byteAryUC:
|     o unsigned char array to hold read in bytes in
|       little endin format
| Output:
|   - Modifies:
|     o byteAryUC to have next bytes (in little endin)
|     o zipMaskUC in fileSTPtr to be 1 (first bit)
|     o zipByteUC to be on byte after byteSI
|     o zipFILE in fileSTPtr to be on byte after zipByteUC
|   - Returns:
|     o 0 for success
|     o def_eofErr_inflate for early EOF
\-------------------------------------------------------*/
signed char
getBytes_inflate(
   struct file_inflate *fileSTPtr,
   signed int bytesSI,
   unsigned char *byteAryUC
){
   if(fileSTPtr->eofBl)
      return def_eofErr_inflate;

   fileSTPtr->zipMaskUC = 8; /*always reseting*/

   /*not best way, but ensures the byte order is correct*/
   while(--bytesSI >= 0)
   { /*Loop: get bytes*/
      if(
         ! fread(
            (char *) &byteAryUC[bytesSI],
            sizeof(unsigned char),
            1,
            fileSTPtr->zipFILE
         )
      ) return def_eofErr_inflate;
   } /*Loop: get bytes*/

   return 0;
} /*getBytes_inflate*/

/*-------------------------------------------------------\
| Fun17: getValue_inflate
|   - gets a value from zlib stream and converts to a
|     value using a huffman tree
| Input:
|   - aryUC:
|     o unsigned char array with value to get
|   - lenBl:
|     o 1: using length tree in fileSTPtr
|     o 0: using distance tree in fileSTPtr
|   - errSCPtr:
|     o signed char pointer to get errors
| Output:
|   - Modifies:
|     o zipPosSL in buffsTPtr to be incurmented for next
|       bits
|     o zipAryUC in fileSTPtr to have new buffer (if
|       needed)
|     o zipFILE in fileSTPtr to be on next part of file
|     o errSCPtr:
|       * 0 if had no errors
|       * def_eofErr_infalte if had early EOF
|   - Returns:
|     o value that code represents
\-------------------------------------------------------*/
unsigned short
getValue_inflate(
   struct file_inflate *fileSTPtr, /*has value to decode*/
   signed char lenBl,              /*1: use length tree*/
   signed char *errSCPtr           /*gets errors*/
){
   struct tree_inflate *iterSTPtr = 0;
   *errSCPtr = 0;

   if(lenBl)
      iterSTPtr = fileSTPtr->lenTreeST;
   else
      iterSTPtr = fileSTPtr->distTreeST;

   while(iterSTPtr->valUS == (unsigned short) -1)
   { /*Loop: get and decode value*/
      if( getBit_inflate(fileSTPtr) )
         iterSTPtr = iterSTPtr->oneST;
      else
         iterSTPtr = iterSTPtr->zeroST;

      if(fileSTPtr->eofBl)
         goto eofErr_fun17;
   } /*Loop: get and decode value*/


   *errSCPtr = 0;
   return iterSTPtr->valUS;

   eofErr_fun17:;
      *errSCPtr = def_eofErr_inflate;
      return 0;
} /*getValue_inflate*/

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
   unsigned short *lenAryUS,      /*length of each code*/
   unsigned short lenUS           /*number of codes*/
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun18 TOC:
   '   - converts a dictionary to a huffman tree
   '   o fun18 sec01:
   '     - variable declarations
   '   o fun18 sec02:
   '     - allocate memory and initialize
   '   o fun18 sec03:
   '     - set up code length and next code arrays
   '   o fun18 sec04:
   '     - build huffman tree
   '   o fun18 sec05:
   '     - return and if needed, clean up
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   struct tree_inflate *retHeapSTPtr = 0;
   struct tree_inflate *nodeSTPtr = 0;

   unsigned short usVal = 0;
   unsigned short usLen = 0;
   unsigned short codeUS = 0;

   #define maxBits_fun18 16
   #define maxAry_fun18 (maxBits_fun18 + 1)
      /*maximum number of bits + 1, so max array length*/
   unsigned short bitCntAryUS[maxAry_fun18];
   unsigned short codeAryUS[maxAry_fun18];

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec02:
   ^   - allocate memory and initialize
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   retHeapSTPtr = malloc(sizeof(struct tree_inflate));
   if(! retHeapSTPtr)
      goto memErr_fun18_sec05;
   init_tree_inflate(retHeapSTPtr);

   for(usVal = 0; usVal < maxAry_fun18; ++usVal)
   { /*Loop: initialize arrays*/
      bitCntAryUS[usVal] = 0;
      codeAryUS[usVal] = 0;
   }  /*Loop: initialize arrays*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec03:
   ^   - set up code length and get codes for each length
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   for(usVal = 0; usVal < lenUS; ++usVal)
     ++bitCntAryUS[ lenAryUS[usVal] ];
        /*build code length array*/

   bitCntAryUS[0] = 0;
   for(usVal = 1; usVal <= maxBits_fun18; ++usVal)
   { /*Loop: get codes for tree*/
      codeUS += bitCntAryUS[usVal - 1];
      codeUS <<= 1;
      codeAryUS[usVal] = codeUS;
   } /*Loop: get codes for tree*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec04:
   ^   - build huffman tree
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   codeUS = 0;

   for(usVal = 0; usVal < lenUS; ++usVal)
   { /*Loop: build huffman tree*/
      if(! lenAryUS[usVal])
         ; /*no length for code, so no code*/

      else
      { /*Else: have length (can add code)*/
         nodeSTPtr = retHeapSTPtr;
         usLen = lenAryUS[usVal];
         codeUS = codeAryUS[usLen];
         --usLen;

         for(; (signed short) usLen >= 0; --usLen)
         { /*Loop: add code to tree*/
            /*usLen will never be greater than signed
            `  short maximum, but codeAryUS can be
            */

            if(codeUS & (1 << usLen))
            { /*If: is a one*/
               if(! nodeSTPtr->oneST)
               { /*If: need to make new node*/
                  nodeSTPtr->oneST =
                     malloc(sizeof(struct tree_inflate));
                  if(! nodeSTPtr->oneST)
                     goto memErr_fun18_sec05;
                  init_tree_inflate(nodeSTPtr->oneST);
               } /*If: need to make new node*/

               nodeSTPtr = nodeSTPtr->oneST;
            } /*If: is a one*/

            else
            { /*Else: is a zero*/
               if(! nodeSTPtr->zeroST)
               { /*If: need to make new node*/
                  nodeSTPtr->zeroST =
                     malloc(sizeof(struct tree_inflate));
                  if(! nodeSTPtr->zeroST)
                     goto memErr_fun18_sec05;
                  init_tree_inflate(nodeSTPtr->zeroST);
               } /*If: need to make new node*/

               nodeSTPtr = nodeSTPtr->zeroST;
            } /*Else: is a zero*/


            if(! usLen)
            { /*If: adding code in*/
               if(nodeSTPtr->valUS != (unsigned short) -1)
                  goto badSymbol_fun18_sec05;
                  /*code already assigned*/
               else
                  nodeSTPtr->valUS = usVal;
            } /*If: adding code in*/
         } /*Loop: add code to tree*/

         ++codeAryUS[lenAryUS[usVal]];
            /*move to lengnths next code*/
      } /*Else: have length (can add code)*/
   } /*Loop: build huffman tree*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun18 Sec05:
   ^   - return and if needed, clean up
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   goto ret_fun18_sec05;

   memErr_fun18_sec05:;
      goto errClean_fun18_sec05;

   badSymbol_fun18_sec05:;
      goto errClean_fun18_sec05;

   errClean_fun18_sec05:;
      if(retHeapSTPtr)
         freeTreeHeap_tree_inflate(retHeapSTPtr);
      retHeapSTPtr = 0;
      goto ret_fun18_sec05;

   ret_fun18_sec05:;
      return retHeapSTPtr;
} /*dictToTree_inflate*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun19 TOC:
   '   - gets huffman trees for inflation
   '   o fun19 sec01:
   '     - variable declarations
   '   o fun19 sec02:
   '     - initialize arrayy
   '   o fun19 sec03:
   '     - get lengths and build length tree for next step
   '   o fun19 sec04:
   '     - build distance tree
   '   o fun19 sec05:
   '     - build length tree
   '   o fun19 sec06:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;

   signed int siLoop = 0;
   signed int repSI = 0;
   unsigned short valUS = 0;

   unsigned short dictLenUS = 0;/*# lenabet and length*/
   unsigned short distUS = 0;   /*number distance codes*/
   unsigned short codeUS = 0;   /*number of code lengths*/

   #define sizeCode_fun19 282
      /*largest value possible is hlit (dictLen)
      `   which is 5 bits + 257 = 16 + 257 = 273
      */
   unsigned short bitLenAryUS[
         sizeCode_fun19 + (sizeof(ulong_ulCp) << 1)
   ];

   #ifdef NOUL
   #else
      /*used to blank arrays (when need to)*/
      ulong_ulCp *bitInitULPtr = 0;
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec02:
   ^   - initialize arrayy
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   #ifdef NOUL
    for(siLoop = 0; siLoop < sizeCode_fun19 + 1; ++siLoop)
       bitLenAryUS[siLoop] = 0;
   #else
      bitInitULPtr = (ulong_ulCp *) bitLenAryUS;
      for(
         siLoop = 0;
         siLoop <
             ((sizeCode_fun19 >> def_shiftULBy_ulCp) << 1)
           + 1;
         ++siLoop
      ) bitInitULPtr[siLoop] = 0;
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec03:
   ^   - get lengths and build length tree for next step
   ^   o fun19 sec03 sub01:
   ^     - get lengths (dictionary, distance, and code)
   ^   o fun19 sec03 sub02:
   ^     - build length/end array for tree
   ^   o fun19 sec03 sub03:
   ^     - build initial length tree
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun19 Sec03 Sub01:
   *   - get lengths (dictionary, distance, and code)
   \*****************************************************/

   dictLenUS =
      getNBits_inflate(
         fileSTPtr,
         5,
         &errSC
      );
   if(errSC)
      goto eofErr_fun19_sec06;
   dictLenUS += 257;

   distUS =
      getNBits_inflate(
         fileSTPtr,
         5,
         &errSC
      );
   if(errSC)
      goto eofErr_fun19_sec06;
   ++distUS;

   codeUS =
      getNBits_inflate(
         fileSTPtr,
         4,
         &errSC
      );
   if(errSC)
      goto eofErr_fun19_sec06;
   codeUS += 4;

   /*****************************************************\
   * Fun19 Sec03 Sub02:
   *   - build length/end array for tree
   \*****************************************************/

   for(siLoop = 0; siLoop < codeUS; ++siLoop)
   { /*Loop: get codes for initial length tree*/
      bitLenAryUS[ codeLenTbl_inflate[siLoop] ] =
         getNBits_inflate(fileSTPtr, 3, &errSC);
      if(errSC)
         goto eofErr_fun19_sec06;
   } /*Loop: get codes for initial length tree*/


   if(fileSTPtr->distTreeST)
      freeTreeHeap_tree_inflate(fileSTPtr->distTreeST);
   fileSTPtr->distTreeST = 0;


   fileSTPtr->distTreeST =
      dictToTree_inflate(
         bitLenAryUS,
         19
      );

   if(! fileSTPtr->distTreeST)
      goto memErr_fun19_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec04:
   ^   - build distance tree
   ^   o fun19 sec04 sub01:
   ^     - initialize bitLenArySS (re-using)
   ^   o fun19 sec04 sub02:
   ^     - get length tree code lengths
   ^   o fun19 sec04 sub03:
   ^     - build length tree
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun19 Sec04 Sub01:
   *   - initialize bitLenArySS (re-using)
   \*****************************************************/

   #ifdef NOUL
    for(siLoop = 0; siLoop < sizeCode_fun19 + 1; ++siLoop)
       bitLenAryUS[siLoop] = 0;
   #else
      bitInitULPtr = (ulong_ulCp *) bitLenAryUS;
      for(
         siLoop = 0;
         siLoop <
             ((sizeCode_fun19 >> def_shiftULBy_ulCp) << 1)
           + 1;
         ++siLoop
      ) bitInitULPtr[siLoop] = 0;
   #endif

   /*****************************************************\
   * Fun19 Sec04 Sub02:
   *   - get length tree code lengths
   \*****************************************************/

   siLoop = 0;
   do{ /*Loop: build length array for real length tree*/
      valUS = getValue_inflate(fileSTPtr, 0, &errSC);
          /*get bit length*/
      if(errSC)
         goto eofErr_fun19_sec06;
      
      if(valUS == (unsigned short) -1)
         goto badSymbol_fun19_sec06;

      else if(valUS <= 15)
      { /*Else If: no repeat of value (symbol)*/
         bitLenAryUS[siLoop++] = valUS;
         continue;
      } /*Else If: no repeat of value (symbol)*/

      else if(valUS == 16)
      { /*Else If: repeating previous code length*/
         repSI = getNBits_inflate( fileSTPtr, 2, &errSC);
         if(errSC)
            goto eofErr_fun19_sec06;
         repSI += 3; /*goes to 3 to 6 times*/
         valUS = bitLenAryUS[siLoop - 1];
      } /*Else If: repeating previous code length*/

      else if(valUS == 17)
      { /*Else If: repeating value (3-10 times)*/
         repSI = getNBits_inflate( fileSTPtr, 3, &errSC);
         if(errSC)
            goto eofErr_fun19_sec06;
         valUS = 0;
         repSI += 3;
      } /*Else If: repeating value (3-10 times)*/

      else if(valUS == 18)
      { /*Else If: repeating value (11 to 138 times)*/
         repSI = getNBits_inflate( fileSTPtr, 7, &errSC);
         if(errSC)
            goto eofErr_fun19_sec06;
         valUS = 0;
         repSI += 11;
      } /*Else If: repeating value (11 to 138 times)*/

      else
         goto badSymbol_fun19_sec06;

      while(repSI--)
         bitLenAryUS[siLoop++] = valUS;
   } while(siLoop < dictLenUS);
    /*Loop: build length array for real length tree*/

   /*****************************************************\
   * Fun19 Sec04 Sub03:
   *   - build length tree
   \*****************************************************/

   if(fileSTPtr->lenTreeST)
      freeTreeHeap_tree_inflate(fileSTPtr->lenTreeST);
   fileSTPtr->lenTreeST = 0;

   fileSTPtr->lenTreeST =
      dictToTree_inflate(
         bitLenAryUS,
         dictLenUS
      );

   if(! fileSTPtr->lenTreeST)
      goto memErr_fun19_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec05:
   ^   - build distance tree
   ^   o fun19 sec05 sub01:
   ^     - initialize bitLenArySS (re-using)
   ^   o fun19 sec05 sub02:
   ^     - get distance tree code lengths
   ^   o fun19 sec05 sub03:
   ^     - build distance tree
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun19 Sec05 Sub01:
   *   - initialize bitLenArySS (re-using)
   \*****************************************************/

   #ifdef NOUL
    for(siLoop = 0; siLoop < sizeCode_fun19 + 1; ++siLoop)
       bitLenAryUS[siLoop] = 0;
   #else
      bitInitULPtr = (ulong_ulCp *) bitLenAryUS;
      for(
         siLoop = 0;
         siLoop <
             ((sizeCode_fun19 >> def_shiftULBy_ulCp) << 1)
           + 1;
         ++siLoop
      ) bitInitULPtr[siLoop] = 0;
   #endif

   /*****************************************************\
   * Fun19 Sec05 Sub02:
   *   - get distance tree code lengths
   \*****************************************************/

   siLoop = 0;
   do{ /*Loop: build distance array for distance tree*/
      valUS = getValue_inflate(fileSTPtr, 0, &errSC);
          /*get bit distance*/
      if(errSC)
         goto eofErr_fun19_sec06;
      
      if(valUS == (unsigned short) -1)
         goto badSymbol_fun19_sec06;

      else if(valUS <= 15)
      { /*Else If: no repeat of value (symbol)*/
         bitLenAryUS[siLoop++] = valUS;
         continue;
      } /*Else If: no repeat of value (symbol)*/

      else if(valUS == 16)
      { /*Else If: repeating previous code distgth*/
         repSI = getNBits_inflate( fileSTPtr, 2, &errSC);
         if(errSC)
            goto eofErr_fun19_sec06;
         repSI += 3; /*goes to 3 to 6 times*/
         valUS = bitLenAryUS[siLoop - 1];
      } /*Else If: repeating previous code distgth*/

      else if(valUS == 17)
      { /*Else If: repeating value (3-10 times)*/
         repSI = getNBits_inflate( fileSTPtr, 3, &errSC);
         if(errSC)
            goto eofErr_fun19_sec06;
         valUS = 0;
         repSI += 3;
      } /*Else If: repeating value (3-10 times)*/

      else if(valUS == 18)
      { /*Else If: repeating value (11 to 138 times)*/
         repSI = getNBits_inflate(fileSTPtr, 7, &errSC);
         if(errSC)
            goto eofErr_fun19_sec06;
         valUS = 0;
         repSI += 11;
      } /*Else If: repeating value (11 to 138 times)*/

      else
         goto badSymbol_fun19_sec06;

      while(repSI--)
         bitLenAryUS[siLoop++] = valUS;
   } while(siLoop < distUS);
    /*Loop: build distance array for distance tree*/

   /*****************************************************\
   * Fun19 Sec05 Sub03:
   *   - build distance tree
   \*****************************************************/

   if(fileSTPtr->distTreeST)
      freeTreeHeap_tree_inflate(fileSTPtr->distTreeST);
   fileSTPtr->distTreeST = 0;

   fileSTPtr->distTreeST =
      dictToTree_inflate(
         bitLenAryUS,
         distUS
      );

   if(! fileSTPtr->distTreeST)
      goto memErr_fun19_sec06;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun19 Sec06:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   distUS = 0;
   goto ret_fun19_sec06;

   memErr_fun19_sec06:;
      distUS = def_memErr_inflate;
      goto ret_fun19_sec06;

   badSymbol_fun19_sec06:;
      distUS = def_badSymbol_inflate;
      goto ret_fun19_sec06;

   eofErr_fun19_sec06:;
      distUS = def_eofErr_inflate;
      goto ret_fun19_sec06;

   ret_fun19_sec06:;
      return (signed char) distUS;
} /*getTree_inflate*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun20 TOC:
   '   - gets huffman fixed trees for fixed inflation
   '   o fun20 sec01:
   '     - variable declarations and initialize
   '   o fun20 sec02:
   '     - build length tree
   '   o fun20 sec03:
   '     - build distance tree
   '   o fun20 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec01:
   ^   - variable declarations and initialize
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed int siLoop = 0;
   signed int indexSI = 0;

   #define sizeLen_fun20 (144 + 112 + 24 + 8)
      /*fixed number of codes*/
   unsigned short aryUS[sizeLen_fun20];

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec02:
   ^   - build length tree
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*this part is untested and comes from
   `   https://pyokagan.name/blog/2019-10-18-zlibinflate/
   `  Hopfully it works
   `  TODO: test this
   */
   for(siLoop = 0; siLoop < 144; ++siLoop)
      aryUS[indexSI++] = 8;
   for(siLoop = 0; siLoop < 112; ++siLoop)
      aryUS[indexSI++] = 9;
   for(siLoop = 0; siLoop < 24; ++siLoop)
      aryUS[indexSI++] = 7;
   for(siLoop = 0; siLoop < 8; ++siLoop)
      aryUS[indexSI++] = 8;

   if(fileSTPtr->lenTreeST)
      freeTreeHeap_tree_inflate(fileSTPtr->lenTreeST);
   fileSTPtr->lenTreeST = 0;

   fileSTPtr->lenTreeST = dictToTree_inflate(aryUS, 286);

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec03:
   ^   - build distance tree
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   indexSI = 0;
   for(siLoop = 0; siLoop < 30; ++siLoop)
      aryUS[indexSI++] = 5;

   if(fileSTPtr->distTreeST)
      freeTreeHeap_tree_inflate(fileSTPtr->distTreeST);
   fileSTPtr->distTreeST = 0;

   fileSTPtr->distTreeST = dictToTree_inflate(aryUS, 30);

   if(! fileSTPtr->distTreeST)
      goto memErr_fun20_sec06;
   
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun20 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   indexSI = 0;
   goto ret_fun20_sec06;

   memErr_fun20_sec06:;
      indexSI = def_memErr_inflate;
      goto ret_fun20_sec06;

   ret_fun20_sec06:;
      return (signed char) indexSI;
} /*getFixedTree_inflate*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun21 TOC:
   '   - inflates data in a block
   '   o fun21 sec01:
   '     - variable declarations
   '   o fun21 sec02:
   '     - check if need initial memory
   '   o fun21 sec03:
   '     - inflate block
   '   o fun21 sec04:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned short valUS = 0;
   signed short repSS = 0;    /*number times to repeat*/
   unsigned short distUS = 0;/*distance to find pattern*/
   unsigned short tmpUS = 0;
   signed long lenCpSL = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec03:
   ^   - inflate block
   ^   o fun21 sec03 sub01:
   ^     - get value & check if need to clear old buffer
   ^   o fun21 sec03 sub02:
   ^     - add single value to internal buffer
   ^   o fun21 sec03 sub03:
   ^     - handle pattern expansion <length,distance>
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun21 Sec03 Sub01:
   *   - get value & check if need to clear old buffer
   \*****************************************************/

   lenCpSL = 0;

   /*not worried about overflow in copy step because
   `  my buffer size is set to the maximum deflate window
   `  size plus two (two is for line breaks) and plus
   `  eight uncounted bytes (for null and ulCp)
   */
   while(
      lenCpSL < def_maxDist_inflate - def_maxLen_inflate
   ){ /*Loop: inflate block*/
      valUS = getValue_inflate(fileSTPtr, 1, errSCPtr);
      if(*errSCPtr)
         goto eofErr_fun21_sec04;

      /**************************************************\
      * Fun21 Sec03 Sub02:
      *   - add single value to internal buffer
      \**************************************************/

      if(valUS <= 255)
      { /*If: literal byte*/

         if(fileSTPtr->buffLenUS >= def_window_inflate)
         { /*If: reached end or array, circle to start*/
            fileSTPtr->buffSizeUS = def_window_inflate;
            fileSTPtr->buffLenUS = 0;
         }  /*If: reached end or array, circle to start*/

         fileSTPtr->buffStr[fileSTPtr->buffLenUS++] = 
            (signed char) valUS;


         #ifdef SIXTEEN_BIT
            if(fileSTPtr->typeSC == def_gz_inflate)
               ;
         #else
            ++fileSTPtr->lenUI;

            if(fileSTPtr->typeSC == def_gz_inflate)
            { /*If: need to get gz crc32 checksum*/
               fileSTPtr->crc32UI =
                  crc32Byte_checkSum(
                     valUS,
                     fileSTPtr->crc32UI
                  );
            } /*If: need to get gz crc32 checksum*/
         #endif
  
         else
         { /*Else: use zlib adler32 checksum*/
            usAdler32Byte_checkSum(
               valUS,
               fileSTPtr->checkSumAryUC
            );
         } /*Else: use zlib adler32 checksum*/


         ++lenCpSL;
      } /*If: literal byte*/

      else if(valUS == 256)
         goto endBlock_fun21_sec04;

      /**************************************************\
      * Fun21 Sec03 Sub03:
      *   - handle pattern expansion <length,distance>
      *   o fun21 sec03 sub03 cat01:
      *     - get length & number bytes (distance) to rep
      *   o fun21 sec03 sub03 cat02:
      *     - copy inflated pattern into buffer
      \**************************************************/
      
      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun21 Sec03 Sub03 Cat01:
      +   - get length and number bytes (distance) to rep
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      else
      { /*Else: is a <length, backwards> pattern*/
         if(valUS < 265)
            repSS = (signed short) (valUS - 254);
         else if(valUS >= 285)
            repSS = 258;
         else
         { /*Else: need to read bits to get length*/
            repSS =
               (signed short)
               getNBits_inflate(
                  fileSTPtr,
                  (valUS - 261) >> 2,
                  errSCPtr
               ); /*get bits in pattern*/
            if(*errSCPtr)
               goto eofErr_fun21_sec04;
            repSS += lenExtraTbl_inflate[valUS - 265];
         } /*Else: need to read bits to get length*/

         distUS =
            getValue_inflate(fileSTPtr, 0, errSCPtr);
         if(*errSCPtr)
            goto eofErr_fun21_sec04;

         if(distUS > 3)
         { /*If: distance is encoded in file*/
            tmpUS =
               getNBits_inflate(
                  fileSTPtr,
                  (distUS - 2) >> 1,
                  errSCPtr
               ); /*get distance to replicate*/
            if(*errSCPtr)
               goto eofErr_fun21_sec04;

				distUS = distExtraTbl_inflate[distUS - 4];
            distUS += tmpUS;
         } /*If: distance is encoded in file*/

         /*++++++++++++++++++++++++++++++++++++++++++++++\
         + Fun21 Sec03 Sub03 Cat02:
         +   - copy inflated pattern into buffer
         \++++++++++++++++++++++++++++++++++++++++++++++*/

         lenCpSL += repSS;

         if(fileSTPtr->buffLenUS > fileSTPtr->buffSizeUS)
            fileSTPtr->buffSizeUS = fileSTPtr->buffLenUS;

         if(fileSTPtr->buffLenUS > distUS)
            distUS = fileSTPtr->buffLenUS - distUS - 1;

         else
         { /*Else: copying starts at opposite end*/
            distUS -= fileSTPtr->buffLenUS;
            distUS = fileSTPtr->buffSizeUS - distUS - 1;
         } /*Else: copying starts at opposite end*/


         if(
              fileSTPtr->buffLenUS + 1
            > fileSTPtr->buffSizeUS - repSS
         ){ /*If: will increase buffer size*/
            if(
                  fileSTPtr->buffLenUS + 1
               >= def_window_inflate - repSS
            ) fileSTPtr->buffSizeUS = def_window_inflate;

            else
               fileSTPtr->buffSizeUS =
                  fileSTPtr->buffLenUS + repSS + 1;
         }  /*If: will increase buffer size*/

         #ifdef SIXTEEN_BIT
         #else
            fileSTPtr->lenUI += repSS;
         #endif

         while(repSS--)
         { /*Loop: copy bytes*/
           if(distUS >= fileSTPtr->buffSizeUS)
               distUS = 0;

           if(fileSTPtr->buffLenUS >= def_window_inflate)
              fileSTPtr->buffLenUS = 0;


            #ifdef SIXTEEN_BIT
               if(fileSTPtr->typeSC == def_gz_inflate)
                  ;
            #else
               if(fileSTPtr->typeSC == def_gz_inflate)
               { /*If: need to get gz crc32 checksum*/
                  fileSTPtr->crc32UI =
                     crc32Byte_checkSum(
                        fileSTPtr->buffStr[distUS],
                        fileSTPtr->crc32UI
                     );
               } /*If: need to get gz crc32 checksum*/
            #endif
  
            else
            { /*Else: use zlib adler32 checksum*/
               usAdler32Byte_checkSum(
                  fileSTPtr->buffStr[distUS],
                  fileSTPtr->checkSumAryUC
               );
            } /*Else: use zlib adler32 checksum*/


            fileSTPtr->buffStr[fileSTPtr->buffLenUS++] =
              fileSTPtr->buffStr[distUS++];
         } /*Loop: copy bytes*/
      } /*Else: is a <length, backwards> pattern*/
   } /*Loop: inflate block*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun21 Sec04:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   *errSCPtr = 0;
   goto ret_fun21_sec04;

   endBlock_fun21_sec04:;
      *errSCPtr = def_endBlock_inflate;
      goto ret_fun21_sec04;

   eofErr_fun21_sec04:;
      *errSCPtr = def_eofErr_inflate;
      goto ret_fun21_sec04;

   ret_fun21_sec04:;
      if(fileSTPtr->buffLenUS >= def_window_inflate)
      { /*If: reached end or array, circle to start*/
            fileSTPtr->buffSizeUS = def_window_inflate;
            fileSTPtr->buffLenUS = 0;
      }  /*If: reached end or array, circle to start*/

      else if(
         fileSTPtr->buffLenUS > fileSTPtr->buffSizeUS
      ) fileSTPtr->buffSizeUS = fileSTPtr->buffLenUS;

      fileSTPtr->buffStr[fileSTPtr->buffSizeUS] = '\0';
         /*always make sure is c-string format*/
      return lenCpSL;
} /*block_inflate*/

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
){
   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec01:
   ^   - variable declarations and find bytes to copy
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   unsigned short tmpUS = 0;
   unsigned short cpUS = 0;
   unsigned char breakAryUC[2];

   if(! bytesUS)
      return 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun22 Sec02:
   ^   - loop copy bytes
   ^   o fun22 sec02 sub01:
   ^     - copy till line ending + start copy loop
   ^   o fun22 sec02 sub02:
   ^     - add in line ending to output buffer
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun22 Sec02 Sub01:
   *   - copy till line ending + start copy loop
   \*****************************************************/

   while(bytesUS)
   { /*Loop: copy bytes from buffer*/
      if(cpUS == def_window_inflate)
         break; /*nothing left to copy*/
      else if(
         fileSTPtr->buffPosUS == fileSTPtr->buffLenUS
      ) break; /*nothing left to copy*/
      else if(
         fileSTPtr->buffPosUS >= fileSTPtr->buffSizeUS
      ) fileSTPtr->buffPosUS = 0;
        /*need to restart at start of array*/

      if(! fileSTPtr->binBl)
      { /*If: working with text data*/
         tmpUS =
             (unsigned short)
             endLine_ulCp(
                &fileSTPtr->buffStr[fileSTPtr->buffPosUS]
             );

         if(tmpUS > bytesUS)
            tmpUS = bytesUS;

         if(
            fileSTPtr->buffLenUS < fileSTPtr->buffPosUS
         ) ; /*else if: need to circle around still*/

         else if(
              fileSTPtr->buffPosUS
            > fileSTPtr->buffLenUS - tmpUS
         ){ /*If: moving past buffer limit*/
            tmpUS = fileSTPtr->buffLenUS;
            tmpUS -= fileSTPtr->buffPosUS;
         }  /*If: moving past buffer limit*/

         if(
               fileSTPtr->buffPosUS
            >= fileSTPtr->buffSizeUS - tmpUS
         ){ /*If: moving past buffer limit*/
            tmpUS = fileSTPtr->buffSizeUS;
            tmpUS -= fileSTPtr->buffPosUS;
         }  /*If: moving past buffer limit*/
      } /*If: working with text data*/

      else
      { /*Else: working with binary data*/
         /*needed because working with a circular array*/
         if(fileSTPtr->buffPosUS <= fileSTPtr->buffLenUS)
            tmpUS = fileSTPtr->buffLenUS;
         else
            tmpUS = fileSTPtr->buffSizeUS;

         tmpUS -= fileSTPtr->buffPosUS;
         tmpUS = min_genMath(bytesUS, tmpUS);

         if(tmpUS > bytesUS)
         { /*If: finished*/
            tmpUS = bytesUS;
            bytesUS = 0;
         } /*If: finished*/
      } /*Else: working with binary data*/

      bytesUS -= tmpUS;

      if(tmpUS)
      { /*If: have bytes to copy (not at line end)*/
         cpLen_ulCp(
            &outStr[posSL],
            &fileSTPtr->buffStr[fileSTPtr->buffPosUS],
            tmpUS
         );
      } /*If: have bytes to copy (not at line end)*/

      posSL += tmpUS;/*position at in buffer*/
      cpUS += tmpUS; /*bytes read in*/
      fileSTPtr->buffPosUS += tmpUS;

      if(fileSTPtr->buffPosUS == fileSTPtr->buffLenUS)
         break;
      else if(! fileSTPtr->buffStr[fileSTPtr->buffPosUS])
         fileSTPtr->buffPosUS = 0; /*end of buffer*/

      /**************************************************\
      * Fun22 Sec02 Sub02:
      *   - add in line ending to output buffer
      \**************************************************/

      if(fileSTPtr->binBl)
         continue;
      else if(bytesUS <= 0)
         break; /*nothing left to copy*/


      if(fileSTPtr->buffLenUS > fileSTPtr->buffPosUS)
         tmpUS =
            fileSTPtr->buffLenUS - fileSTPtr->buffPosUS;
      else
      { /*Else: need to wrap arround array*/
         tmpUS = fileSTPtr->buffSizeUS;
         tmpUS -= fileSTPtr->buffPosUS;
         tmpUS += fileSTPtr->buffLenUS;
      } /*Else: need to wrap arround array*/


      if( ! (newLineBl & 4) && bytesUS < 2 )
         break; /*not enough room for line break*/
      else if(tmpUS > 1)
         ; /*at worst case; have full line break*/
      else if(! fileSTPtr->finalBl)
         break; /*not last block and not full break*/
      else if(
         fileSTPtr->blockFlagSC != def_noBlock_inflate
      ) break; /*have not finished reading final block*/
      /*Else have enough bytes in buffer and either
      `  2 bytes in tmpUS or am at end of final
      `  block (no more bytes)
      */

      /*find size of line break*/
      if(str_endLine[1])
         tmpUS = 2;
      else
         tmpUS = 1;


      if(tmpUS <= bytesUS)
      { /*Else If: have room for line ending*/
         tmpUS = fileSTPtr->buffPosUS;
         breakAryUC[0] = fileSTPtr->buffStr[tmpUS];

         if(tmpUS + 1 < fileSTPtr->buffSizeUS)
            breakAryUC[1] = fileSTPtr->buffStr[tmpUS];
         else if(
            fileSTPtr->buffLenUS < fileSTPtr->buffSizeUS
         ) breakAryUC[1] = fileSTPtr->buffStr[0];
         else
            breakAryUC[1] = fileSTPtr->buffStr[tmpUS];
         

         if(
            breakAryUC[0] == '\n' && breakAryUC[1] == '\r'
         ) tmpUS = 2;
         else if(
            breakAryUC[0] == '\r' && breakAryUC[1] == '\n'
         ) tmpUS = 2;
         else if(breakAryUC[0] == '\n')
            tmpUS = 1;
         else if(breakAryUC[0] == '\r')
            tmpUS = 1;
         else
            tmpUS = 0;

         fileSTPtr->buffPosUS += tmpUS;
         if(fileSTPtr->buffPosUS == fileSTPtr->buffLenUS)
            ; /*at end of buffer*/
         else if(
            fileSTPtr->buffPosUS >= fileSTPtr->buffSizeUS
         ) fileSTPtr->buffPosUS -= fileSTPtr->buffSizeUS;


         if(tmpUS)
         { /*If: had line break*/
            outStr[posSL++] = str_endLine[0];
            ++cpUS;
            --bytesUS;

            if(str_endLine[1] != '\0')
            { /*If: have a two character line ending*/
               outStr[posSL++] = str_endLine[1];
               ++cpUS;
               --bytesUS;
            } /*If: have a two character line ending*/

            if(newLineBl)
               break;
         } /*If: had line break*/
      } /*Else If: have room for line ending*/

      else
         break; /*not enough room for line ending*/
   } /*Loop: copy bytes from buffer*/

   outStr[posSL] = '\0';
   return cpUS;
} /*cpBytes_inflate*/

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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun23 TOC:
   '   - inflate blocks until have required input
   '   o fun23 sec01:
   '     - variable declarations
   '   o fun23 sec02:
   '     - copy bytes remaining from last read
   '   o fun23 sec03:
   '     - get bytes
   '   o fun23 sec04:
   '     - return
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char newBlockBl = 0; /*1:starting on new block*/
   unsigned short lenUS = 0;
   unsigned short compLenUS = 0;
   signed char lineBl = (endBl & 1);
   signed char extraBytesBl = endBl & 4;

   signed long cpLenSL = 0;
   unsigned short tmpUS = 0;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec02:
   ^   - have bytes in last bock to copy over
   ^   o fun23 sec02 sub01:
   ^     - if user wanted not bytes or need more bytes
   ^   o fun23 sec02 sub02:
   ^     - copy uncopied bytes to buffer
   ^   o fun23 sec02 sub03:
   ^     - if no copied bytes, make sure not end of file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun23 Sec02 Sub01:
   *   - if user wanted not bytes or need more bytes
   \*****************************************************/

   if( bytesSL <= 0)
   { /*If: not bytes to copy*/
      outStr[0] = '\0'; /*user wanted no bytes*/
      if(! (endBl & 8) )
         return 0; /*not refilling on edge case*/
   } /*If: not bytes to copy*/


   if(
         fileSTPtr->buffLenUS - fileSTPtr->buffPosUS == 1
      && ! fileSTPtr->finalBl
   ) ; /*likely at new line, fill buffer*/

   /*****************************************************\
   * Fun23 Sec02 Sub02:
   *   - copy uncopied bytes to buffer
   \*****************************************************/

   else if(endBl & 8)
      ;
   else if(fileSTPtr->buffPosUS != fileSTPtr->buffLenUS)
   { /*Else If: have bytes to copy*/
      copyBytes_fun23_sec02_sub01:;

      if(bytesSL > (unsigned short) -1)
         tmpUS = ((unsigned short) -1) - 2;
      else
         tmpUS = (unsigned short) bytesSL;

      tmpUS =
         cpBytes_inflate(
            fileSTPtr,
            outStr,
            cpLenSL,
            lineBl | extraBytesBl,
            tmpUS
         ); /*copy any extra bytes in buffer*/


      cpLenSL += tmpUS;
      bytesSL -= tmpUS;

      if(bytesSL <= 0)
      { /*If: read in all requested bytes*/
         if(! fileSTPtr->finalBl)
            goto noErr_fun23_sec04; /*not on final block*/

         else if(
            fileSTPtr->buffPosUS != fileSTPtr->buffLenUS
         ) goto noErr_fun23_sec04; /*not at end of block*/

         else if(
            fileSTPtr->blockFlagSC != def_noBlock_inflate
         ) goto noErr_fun23_sec04; /*more block to read*/

         else
            goto eof_fun23_sec04; /*finshed final block*/
      } /*If: read in all requested bytes*/

      else if(
         fileSTPtr->buffPosUS == fileSTPtr->buffLenUS
      ){ /*Else If: at end of read in buffer*/
         if(lineBl && outStr[cpLenSL - 1] == '\n')
            goto noErr_fun23_sec04; /*hit end of line*/
         else if(lineBl && outStr[cpLenSL - 1] == '\r')
            goto noErr_fun23_sec04; /*hit end of line*/
         else if( bytesSL == 1 && ! (endBl & 4) )
            goto noErr_fun23_sec04; /*buffer is full*/
         else if(
            fileSTPtr->blockFlagSC != def_noBlock_inflate
         ) ; /*not at end of block*/
         else if(fileSTPtr->finalBl)
            goto eof_fun23_sec04; /*finshed final block*/
         else if(endBl & 2)
            goto noErr_fun23_sec04; /*got block*/
         /*else; grab more bytes*/
      }  /*Else If: at end of read in buffer*/

      else if(lineBl && outStr[cpLenSL - 1] == '\n')
         goto noErr_fun23_sec04; /*hit end of line*/
      else if(lineBl && outStr[cpLenSL - 1] == '\r')
         goto noErr_fun23_sec04; /*hit end of line*/
      else if(! bytesSL)
         goto noErr_fun23_sec04; /*buffer is full*/
      else if(bytesSL >= 2)
         ; /*grab more bytes*/
      else if(endBl & 4)
         ;
      else
         goto noErr_fun23_sec04;
         /*not enough room to ensure got line ending*/
      /*else; grab more bytes*/
   } /*Else If: have bytes to copy*/

   /*****************************************************\
   * Fun23 Sec02 Sub03:
   *   - if no copied bytes, make sure not at end of file
   \*****************************************************/

   else if(! fileSTPtr->finalBl)
      ; /*not on final block*/
   else if(fileSTPtr->blockFlagSC == def_noBlock_inflate)
      goto eof_fun23_sec04; /*finshed final block*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec03:
   ^   - get bytes
   ^   o fun23 sec03 sub01:
   ^     - get block type (if unkown) + start block loop
   ^   o fun23 sec03 sub02:
   ^     - no compression blocks read in
   ^   o fun23 sec03 sub03:
   ^     - get tree for compressed block
   ^   o fun23 sec03 sub04:
   ^     - inflate compressed block
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun23 Sec03 Sub01:
   *   - get block type (if unkown) + start block loop
   \*****************************************************/

   *errSCPtr = 0;

   if(fileSTPtr->blockFlagSC == def_noBlock_inflate)
   { /*If: first time or a new block*/
      newBlockBl = 1;
      fileSTPtr->bytesLeftUS = 0;
         /*restart uncompressed byte block count*/

      fileSTPtr->finalBl =
         (signed char) getBit_inflate(fileSTPtr);
      if(fileSTPtr->eofBl)
      { /*If: hit EOF*/
         if(fileSTPtr->finalBl)
            goto eof_fun23_sec04;
         else
            goto eofErr_fun23_sec04;
      } /*If: hit EOF*/


      fileSTPtr->blockFlagSC =
         (signed char) getBit_inflate(fileSTPtr);
      if(fileSTPtr->eofBl)
         goto eofErr_fun23_sec04;

      fileSTPtr->blockFlagSC |=
         ( getBit_inflate(fileSTPtr) << 1 );
      if(fileSTPtr->eofBl)
         goto eofErr_fun23_sec04;
   } /*If: first time or a new block*/

   /**************************************************\
   * Fun23 Sec03 Sub02:
   *   - read in no compression blocks
   \**************************************************/

   if(! fileSTPtr->blockFlagSC)
   { /*If: no compression*/
      if(newBlockBl)
      { /*If: at start of block*/
         newBlockBl = 0;

         if(
            getBytes_inflate(
               fileSTPtr,
               2,
               (unsigned char *) &lenUS
            )
         ) goto eofErr_fun23_sec04;

         if(
            getBytes_inflate(
               fileSTPtr,
               2,
               (unsigned char *) &compLenUS
            )
         ) goto eofErr_fun23_sec04;

         fileSTPtr->bytesLeftUS = lenUS;
      } /*If: at start of block*/


      tmpUS =
         min_genMath(
            fileSTPtr->bytesLeftUS,
            def_window_inflate
         );

      fileSTPtr->buffLenUS =
          (unsigned short)
          fread(
            (char *) fileSTPtr->buffStr,
            sizeof(signed char),
            tmpUS,
            fileSTPtr->zipFILE
         );

      fileSTPtr->zipMaskUC = 8;
      fileSTPtr->bytesLeftUS -= fileSTPtr->buffLenUS;

      if(fileSTPtr->buffLenUS < tmpUS)
         goto eofErr_fun23_sec04;
      else if(! fileSTPtr->bytesLeftUS)
         fileSTPtr->blockFlagSC = def_noBlock_inflate;

      fileSTPtr->buffPosUS = 0; /*reset for copying*/
      goto copyBytes_fun23_sec02_sub01;
   } /*If: no compression*/

   /**************************************************\
   * Fun23 Sec04 Sub03:
   *   - get tree for compressed block
   \**************************************************/

   else if(newBlockBl)
   { /*Else If: getting a new block*/
      newBlockBl = 0;

      if(fileSTPtr->blockFlagSC == 1)
      { /*If: fixed tree block*/
         if( getFixedTree_inflate(fileSTPtr) )
            goto memErr_fun23_sec04;
      } /*If: fixed tree block*/


      else if(fileSTPtr->blockFlagSC == 2)
      { /*Else If: dynamic block*/
         *errSCPtr = getTree_inflate(fileSTPtr);

         if(! *errSCPtr)
            ;
         else if(*errSCPtr == def_memErr_inflate)
            goto memErr_fun23_sec04;
         else if(*errSCPtr == def_badSymbol_inflate)
            goto badSymbol_fun23_sec04;
         else if(*errSCPtr == def_eofErr_inflate)
            goto eofErr_fun23_sec04;
      } /*Else If: dynamic block*/

      else
         goto badBlock_fun23_sec04;
   } /*Else If: getting a new block*/

   /**************************************************\
   * Fun23 Sec05 Sub04:
   *   - inflate compressed block
   \**************************************************/

   tmpUS = block_inflate(fileSTPtr, errSCPtr);

   if(! *errSCPtr)
      ;
   else if(*errSCPtr == def_endBlock_inflate)
      fileSTPtr->blockFlagSC = def_noBlock_inflate;
      /*so know at end of block later*/
   else if(*errSCPtr == def_memErr_inflate)
      goto memErr_fun23_sec04;
   else if(*errSCPtr == def_eofErr_inflate)
      goto eofErr_fun23_sec04;

   goto copyBytes_fun23_sec02_sub01;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun23 Sec04:
   ^   - return
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   noErr_fun23_sec04:;
      *errSCPtr = 0;
      goto ret_fun23_sec04;

   eof_fun23_sec04:;
      *errSCPtr = def_eof_inflate;
      goto ret_fun23_sec04;

   memErr_fun23_sec04:;
      *errSCPtr = def_memErr_inflate;
      cpLenSL = -1;
      goto ret_fun23_sec04;

   badSymbol_fun23_sec04:;
      *errSCPtr = def_badSymbol_inflate;
      cpLenSL = -1;
      goto ret_fun23_sec04;

   badBlock_fun23_sec04:;
      *errSCPtr = def_badBlock_inflate;
      cpLenSL = -1;
      goto ret_fun23_sec04;

   eofErr_fun23_sec04:;
      *errSCPtr = def_eofErr_inflate;
      cpLenSL = -1;
      goto ret_fun23_sec04;

   ret_fun23_sec04:;
      return cpLenSL;
} /*get_inflate*/

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
){
   signed char errSC = 0;

   if(fileSTPtr->buffPosUS > fileSTPtr->buffLenUS) ;

   else if(
      fileSTPtr->buffPosUS > fileSTPtr->buffLenUS - 1
   ){ /*Else If: need to grab more bytes*/
      get_inflate(fileSTPtr, 0, 8, getArySC, &errSC);

      if(errSC)
         return errSC * -1;
   }  /*Else If: need to grab more bytes*/

   if(fileSTPtr->binBl)
   { /*If: reading binary data*/
      errSC = 1;
      getArySC[0] =
         fileSTPtr->buffStr[fileSTPtr->buffPosUS];

      if(fileSTPtr->buffPosUS >= fileSTPtr->buffSizeUS)
         fileSTPtr->buffPosUS = 0;
      else
         ++fileSTPtr->buffPosUS;
   } /*If: reading binary data*/

   else if(
      fileSTPtr->buffStr[fileSTPtr->buffPosUS] == '\n'
   ){ /*Else If: have line break*/
      if(fileSTPtr->buffPosUS >= fileSTPtr->buffSizeUS)
         fileSTPtr->buffPosUS = 0;
      else
         ++fileSTPtr->buffPosUS;

      if(fileSTPtr->buffStr[fileSTPtr->buffPosUS] == '\r')
      { /*If: have a two line break ending*/
         if(fileSTPtr->buffPosUS >= fileSTPtr->buffSizeUS)
            fileSTPtr->buffPosUS = 0;
         else
            ++fileSTPtr->buffPosUS;
      } /*If: have a two line break ending*/

      getArySC[0] = str_endLine[0];
      ++errSC;

      if(str_endLine[1])
      { /*If: using 2 character new line line break*/
         ++errSC;
         getArySC[1] = str_endLine[1];
         getArySC[2] = 0;
      } /*If: using 2 character new line line break*/

      else
         getArySC[1] = 0;
   } /*Else If: have line break*/

   else if(
      fileSTPtr->buffStr[fileSTPtr->buffPosUS] == '\r'
   ){ /*Else If: have carriage return line break*/
      if(fileSTPtr->buffPosUS >= fileSTPtr->buffSizeUS)
         fileSTPtr->buffPosUS = 0;
      else
         ++fileSTPtr->buffPosUS;

      if(fileSTPtr->buffStr[fileSTPtr->buffPosUS] == '\n')
      { /*If: have a two line break ending*/
         if(fileSTPtr->buffPosUS >= fileSTPtr->buffSizeUS)
            fileSTPtr->buffPosUS = 0;
         else
            ++fileSTPtr->buffPosUS;
      } /*If: have a two line break ending*/

      getArySC[0] = str_endLine[0];
      ++errSC;

      if(str_endLine[1])
      { /*If: using 2 character line break*/
         ++errSC;
         getArySC[1] = str_endLine[1];
         getArySC[2] = 0;
      } /*If: using 2 character line break*/

      else
         getArySC[1] = 0;
   }  /*Else If: have carriage return line break*/

   else
   { /*Else: non-line break character*/
      errSC = 1;
      getArySC[0] =
         fileSTPtr->buffStr[fileSTPtr->buffPosUS];
      getArySC[1] = 0;

      if(fileSTPtr->buffPosUS >= fileSTPtr->buffSizeUS)
         fileSTPtr->buffPosUS = 0;
      else
         ++fileSTPtr->buffPosUS;
   } /*Else: non-line break character*/

   return errSC;
} /*getc_file_inflate*/

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
){
   if(! fileSTPtr->buffPosUS)
      fileSTPtr->buffPosUS = fileSTPtr->buffSizeUS - 1;
   --fileSTPtr->buffPosUS;
   fileSTPtr->buffStr[fileSTPtr->buffPosUS] = ungetSC;
} /*unget_file_inflate*/

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
|     o def_eof_inflate if at end of file
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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun26 TOC:
   '   - move to next deflate/gz entry in the file
   '   o fun26 sec01:
   '     - variable declarations
   '   o fun26 sec02:
   '     - get past footers for deflate and gz
   '   o fun26 sec03:
   '     - get headers
   '   o fun26 sec04:
   '     - get and check CINo fO block (byte) & checksum
   '   o fun26 sec05:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec01:
   ^   - variable declarations
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC = 0;

   unsigned char cmfUC = 0;    /*is 1st byte*/
   unsigned char cinfoUC = 0;  /*is 2nd byte*/
   unsigned char fileFlagsUC = 0; /*gzip*/
   signed short winLenSS = 0;
   unsigned short lenUS = 0;
      /*length of uncompressed block*/

   unsigned char tmpAryUC[8];

   #ifdef SIXTEEN_BIT
   #else
      unsigned int crc32UI = 0;
      unsigned int fileLenUI = 0;
   #endif

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec02:
   ^   - get past footers for deflate and gz
   ^   o fun26 sec02 sub01:
   ^     - gz format; get ending head
   ^   o fun26 sec02 sub02:
   ^     - zlib format; get ending head
   ^   o fun26 sec02 sub03:
   ^     - all formats; check if at end of file
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun26 Sec02 Sub01:
   *   - gz format; get ending head
   \*****************************************************/

   if(! newBl)
   { /*If: need to clear last entry*/

      if(fileSTPtr->typeSC == def_gz_inflate)
      { /*If: this was a gzip file*/
         if(
            ! fread(
               (char *) &crc32UI,
               sizeof(int),
               1,
               fileSTPtr->zipFILE
            )
         ) goto eofErr_fun26_sec05;

         if(
            ! fread(
               (char *) &fileLenUI,
               sizeof(int),
               1,
               fileSTPtr->zipFILE
            )
         ) goto eofErr_fun26_sec05;


         fileSTPtr->crc32UI =
           crc32Finish_checkSum(fileSTPtr->crc32UI);
         if(crc32UI != fileSTPtr->crc32UI)
            goto badCheck_fun26_sec05;
         if(fileLenUI != fileSTPtr->lenUI)
            goto badCheck_fun26_sec05;


         blank_file_inflate(fileSTPtr, 0);
         fileSTPtr->typeSC = def_gz_inflate;
      } /*If: this was a gzip file*/

      /**************************************************\
      * Fun26 Sec02 Sub02:
      *   - zlib format; get ending head
      \**************************************************/

      else
      { /*Else: zlib format*/
         usAdler32Finish_checkSum(
            fileSTPtr->checkSumAryUC
         ); /*addler needs shifting at end*/

         /*need to get zlib adler32 checksum at end*/
         if(firstBytesArySC)
         { /*If: user input first two bytes*/
            tmpAryUC[3] =
               (unsigned char) firstBytesArySC[0];
            tmpAryUC[2] =
               (unsigned char) firstBytesArySC[1];
            firstBytesArySC = 0;

            errSC =
               getBytes_inflate(
                  fileSTPtr,
                  2,
                  tmpAryUC
               ); /*adler 32 checksum for deflate*/
         } /*If: user input first two bytes*/

         else
            errSC =
               getBytes_inflate(
                  fileSTPtr,
                  4,
                  tmpAryUC
               ); /*adler 32 checksum for deflate*/

         if(errSC == def_eofErr_inflate)
            goto eofErr_fun26_sec05;
         else if(
            tmpAryUC[0] != fileSTPtr->checkSumAryUC[0]
         ) goto badCheck_fun26_sec05;
         else if(
            tmpAryUC[1] != fileSTPtr->checkSumAryUC[1]
         ) goto badCheck_fun26_sec05;
         else if(
            tmpAryUC[2] != fileSTPtr->checkSumAryUC[2]
         ) goto badCheck_fun26_sec05;
         else if(
            tmpAryUC[3] != fileSTPtr->checkSumAryUC[3]
         ) goto badCheck_fun26_sec05;

         blank_file_inflate(fileSTPtr, 0);
         fileSTPtr->typeSC = def_zlib_inflate;
      } /*Else: zlib format*/

      /**************************************************\
      * Fun26 Sec02 Sub03:
      *   - all formats; check if at end of file
      \**************************************************/

      errSC = getBytes_inflate(fileSTPtr, 2, tmpAryUC);
         /*get gz header magic number*/

      if(errSC == def_eofErr_inflate)
         goto eof_fun26_sec05;
      else
         fileSTPtr->zipMaskUC = 8;
   } /*If: need to clear last entry*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec03:
   ^   - get headers
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*****************************************************\
   * Fun26 Sec03 Sub01:
   *   - get gz header (if relavent)
   *   o fun26 sec03 sub01 cat01:
   *     - read in required part of header
   *   o fun26 sec03 sub01 cat02:
   *     - get file data type (binary or text)
   *   o fun26 sec03 sub01 cat03:
   *     - skip extra entry (if present)
   *   o fun26 sec03 sub01 cat04:
   *     - get file name (if present)
   *   o fun26 sec03 sub01 cat05:
   *     - skip comment (if present)
   *   o fun26 sec03 sub01 cat06:
   *     - skip crc16 (if present)
   \*****************************************************/

   /*++++++++++++++++++++++++++++++++++++++++++++++++++++\
   + Fun26 Sec03 Sub01 Cat01:
   +   - read in required part of header
   \++++++++++++++++++++++++++++++++++++++++++++++++++++*/

   if(fileSTPtr->typeSC == def_gz_inflate)
   { /*If: gz file type*/
      if(firstBytesArySC)
      { /*If: already read in the first two bytes*/
         tmpAryUC[0] = firstBytesArySC[1];
         tmpAryUC[1] = firstBytesArySC[0];
         firstBytesArySC = 0;
      } /*If: already read in the first two bytes*/

      else if(! newBl)
         ; /*already read in bytes in footer read step*/
      else
         errSC = getBytes_inflate(fileSTPtr, 2, tmpAryUC);
         /*get gz header magic number*/

      if(errSC == def_eofErr_inflate)
         goto eofErr_fun26_sec05;
      else if(tmpAryUC[1] != 0x1f)
         goto badType_fun26_sec05;
      else if(tmpAryUC[0] != 0x8b)
         goto badType_fun26_sec05;
         /*backwards to account for little endin*/
 

      errSC = getBytes_inflate(fileSTPtr, 1, &cmfUC);
         /*get compression method*/
      if(errSC == def_eofErr_inflate)
         goto eofErr_fun26_sec05;
      else if(cmfUC != 8)
         goto badType_fun26_sec05;


      errSC =
         getBytes_inflate(fileSTPtr, 1, &fileFlagsUC);
         /*get file flags*/
      if(errSC == def_eofErr_inflate)
         goto eofErr_fun26_sec05;


      if(
         ! fread(
            (char *) fileSTPtr->timeAryUC,
            sizeof(unsigned char),
            4,
            fileSTPtr->zipFILE
         ) /*get time stamp*/
      ) goto eofErr_fun26_sec05;

      fileSTPtr->zipMaskUC = 8;

      errSC = getBytes_inflate(fileSTPtr, 1, &cmfUC);
         /*get compression flags (ignore)*/
      if(errSC == def_eofErr_inflate)
         goto eofErr_fun26_sec05;

      errSC = getBytes_inflate( fileSTPtr, 1, &cmfUC);
         /*get operating system (ignore)*/
      if(errSC == def_eofErr_inflate)
         goto eofErr_fun26_sec05;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun26 Sec03 Sub01 Cat02:
      +   - get file data type (binary or text)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*1 is text based, I assume this, which is a
      `  mistake
      */
      if(fileFlagsUC & 1)
         fileSTPtr->binBl = 1;
      else
         fileSTPtr->binBl = 0;

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun26 Sec03 Sub01 Cat03:
      +   - skip extra entry (if present)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(fileFlagsUC & 4)
      { /*If: have extra entry*/
         errSC =
            getBytes_inflate(
               fileSTPtr,
               2,
              (unsigned char *) &lenUS
            ); /*get the length of the subfeilds*/
         if(errSC == def_eofErr_inflate)
            goto eofErr_fun26_sec05;

         fileSTPtr->zipMaskUC = 8;
         lenUS = usToBig_endin(lenUS);
            /*on bigendin systems flips to bigendin, but
            `  does nothing on little endin
            */
         if(lenUS <= def_window_inflate)
         { /*If: can read extra in one go*/
            if(
               ! fread(
                  (char *) fileSTPtr->buffStr,
                  sizeof(signed char),
                  lenUS,
                  (FILE *) fileSTPtr->zipFILE
               )
            ) goto eofErr_fun26_sec05;
         } /*If: can read extra in one go*/

         else
         { /*Else: need two reads to get past extra*/
            if(
               ! fread(
                  (char *) fileSTPtr->buffStr,
                  sizeof(signed char),
                  lenUS - def_window_inflate,
                  (FILE *) fileSTPtr->zipFILE
               )
            ) goto eofErr_fun26_sec05;

            if(
               ! fread(
                  (char *) fileSTPtr->buffStr,
                  sizeof(signed char),
                  def_window_inflate,
                  (FILE *) fileSTPtr->zipFILE
               )
            ) goto eofErr_fun26_sec05;
         } /*Else: need two reads to get past extra*/
      } /*If: have extra entry*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun26 Sec03 Sub01 Cat04:
      +   - get file name (if present)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*this method is not the most efficent, but this
      `   step only happens once per file
      */
      if(fileFlagsUC & 8)
      { /*If: have file name in header*/
         fileSTPtr->nameLenUC = (unsigned char) -1;

         do{ /*Loop: get file name*/
            ++fileSTPtr->nameLenUC;
            /*using first rounds overlfow to set to 0*/

            if(
               ! fread(
                  (char *) &fileSTPtr->nameStr[
                              fileSTPtr->nameLenUC
                  ],
                  sizeof(signed char),
                  1,
                  fileSTPtr->zipFILE
               )
            ) goto eofErr_fun26_sec05;
         }while(fileSTPtr->nameStr[fileSTPtr->nameLenUC]);
           /*Loop: get file name*/

         fileSTPtr->zipMaskUC = 8;
      } /*If: have file name in header*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun26 Sec03 Sub01 Cat05:
      +   - skip comment (if present)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      /*this method is not the most efficent, but this
      `   step only happens once per file
      */
      if(fileFlagsUC & 16)
      { /*If: have comment in header*/
         do{ /*Loop: get file name*/
            if(
               ! fread(
                  (char *) &fileSTPtr->zipByteUC,
                  sizeof(unsigned char),
                  1,
                  fileSTPtr->zipFILE
               )
            ) goto eofErr_fun26_sec05;
         } while(fileSTPtr->zipByteUC);

         fileSTPtr->zipMaskUC = 8;
      } /*If: have comment in header*/

      /*+++++++++++++++++++++++++++++++++++++++++++++++++\
      + Fun26 Sec03 Sub01 Cat06:
      +   - skip crc16 (if present)
      \+++++++++++++++++++++++++++++++++++++++++++++++++*/

      if(fileFlagsUC & 2)
      { /*If: crc-16 in header*/
         errSC = getBytes_inflate(fileSTPtr, 2, tmpAryUC);
            /*move past crc16 in header*/
         if(errSC == def_eofErr_inflate)
            goto eofErr_fun26_sec05;
      } /*If: crc-16 in header*/

      /*bits 5 to 7 (32 to 128) are not used (reserved)*/
   } /*If: gz file type*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec03:
   ^   - get and check CMF block (byte)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   else if(fileSTPtr->typeSC == def_zlib_inflate)
   { /*Else If: zlib file*/
      if(firstBytesArySC)
         cmfUC = firstBytesArySC[0];
      else if(! newBl)
         cmfUC = tmpAryUC[0];
      else
      { /*Else: need to read in CMF block*/
         errSC =
            getBytes_inflate(
               fileSTPtr,
               sizeof(unsigned char),
               &cmfUC
            );
         if(errSC == def_eofErr_inflate)
            goto eofErr_fun26_sec05;
      } /*Else: need to read in CMF block*/

      if( (cmfUC & 0xf) != 8 )
         goto badType_fun26_sec05;

      winLenSS = ( (signed short) (cmfUC & 0xf0) ) >> 4;

      if(winLenSS > 7)
         goto badWindow_fun26_sec05;
      /*winLenSS = (signed short) 1 << (winLenSS + 8);*/
         /*do not use window length, assume max is used*/

      /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
      ^ Fun26 Sec04:
      ^   - get and check CINFO block (byte) & checksum
      \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

      if(firstBytesArySC)
      { /*If: user provided first two bytes of header*/
         cinfoUC = firstBytesArySC[1];
         firstBytesArySC = 0;
      } /*If: user provided first two bytes of header*/

      else if(! newBl)
         cinfoUC = tmpAryUC[1];

      else
      { /*Else: need to read in CINFO block*/
         errSC =
            getBytes_inflate(
               fileSTPtr,
               sizeof(unsigned char),
               &cinfoUC
            );
         if(errSC == def_eofErr_inflate)
            goto eofErr_fun26_sec05;
      } /*Else: need to read in CINFO block*/

      if( ((cmfUC << 8) + cinfoUC) % 31 )
         goto badCheck_fun26_sec05;
         /*cmfUC << 8 = cmfUC * 256*/

      if(cinfoUC & 0x10)
         goto dictErr_fun26_sec05;
         /*dictionarys are not supported*/
   } /*Else If: zlib file*/

   else
      goto badType_fun26_sec05; /*no idea file type*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun26 Sec05:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   errSC = 0;
   goto ret_fun26_sec05;

   eof_fun26_sec05:;
      errSC = def_eof_inflate;
      goto ret_fun26_sec05;

   badType_fun26_sec05:;
      errSC = def_badType_inflate;
      goto ret_fun26_sec05;

   badWindow_fun26_sec05:;
      errSC = def_badWindow_inflate;
      goto ret_fun26_sec05;

   badCheck_fun26_sec05:;
      errSC = def_badCheck_inflate;
      goto ret_fun26_sec05;

   dictErr_fun26_sec05:;
      errSC = def_dictErr_inflate;
      goto ret_fun26_sec05;

   eofErr_fun26_sec05:;
      errSC = def_eofErr_inflate;
      goto ret_fun26_sec05;

   ret_fun26_sec05:;
      return errSC;
} /*next_file_inflate*/

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
|       * you should set zipFILE to 0/null. Let
|         freeStack_file_inflate or freeHeap_file_inflate
|         handle closing this.
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
){ /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
   ' Fun27 TOC:
   '   - decompresses zip data
   '   o fun27 sec01:
   '     - variable declarations and add file in
   '   o fun27 Sec02:
   '     - get first two bytes (identifies if gz file)
   '   o fun27 sec03:
   '     - return result
   \~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun27 Sec01:
   ^   - variable declarations and add file in
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   signed char errSC;
   unsigned char tmpAryUC[8];

   if(! zipFILE)
      goto fileErr_fun27_sec03;

   blank_file_inflate(fileSTPtr, 2);/*will close zipFILE*/
   fileSTPtr->zipFILE = zipFILE;

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun27 Sec02:
   ^   - get first two bytes (identifies if gz file)
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   if(firstBytesArySC)
   { /*If: user supplied first two bytes*/
      tmpAryUC[1] = firstBytesArySC[0];
      tmpAryUC[0] = firstBytesArySC[1];
   } /*If: user supplied first two bytes*/

   else
   { /*Else: need to get first two bytes from header*/
      errSC = getBytes_inflate(fileSTPtr, 2, tmpAryUC);
      if(errSC == def_eofErr_inflate)
         goto eofErr_fun27_sec03;
   } /*If: gz file type*/

   if(tmpAryUC[1] == 0x1f && tmpAryUC[0] == 0x8b)
      fileSTPtr->typeSC = def_gz_inflate;
   else
      fileSTPtr->typeSC = def_zlib_inflate;

   /*account for next function flipping these entries
   `  arround
   */
   tmpAryUC[0] ^= tmpAryUC[1];
   tmpAryUC[1] ^= tmpAryUC[0];
   tmpAryUC[0] ^= tmpAryUC[1];

   errSC =
      next_file_inflate(
         fileSTPtr,
         (signed char *) tmpAryUC,
         1 /*first read in*/
      ); /*read in rest of header*/

   /*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
   ^ Fun27 Sec03:
   ^   - return result
   \<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

   /*error type already set*/
   goto ret_fun27_sec03;

   fileErr_fun27_sec03:;
      errSC = def_fileErr_inflate;
      goto ret_fun27_sec03;

   eofErr_fun27_sec03:;
      errSC = def_eofErr_inflate;
      goto ret_fun27_sec03;

   ret_fun27_sec03:;
      return errSC;
} /*add_file_inflate*/

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
