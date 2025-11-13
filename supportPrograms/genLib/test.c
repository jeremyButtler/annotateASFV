#include <stdio.h>
#include "mkPng.h"
#include "endLine.h"

int
main(
   void
){
   signed int ySI = 0;
   signed char errSC = 0;
   struct st_mkPng *pngSTPtr = mk_st_mkPng(0, 0, 0);
   FILE *outFILE = fopen("del.png", "wb");

   if(! pngSTPtr)
      goto memErr_main;
   if(! outFILE)
      goto fileErr_main;

   for(ySI = 0; ySI < 100; ++ySI)
   { /*Loop: pixel by pixel add*/
       addPixel_st_mkPng(pngSTPtr, 14, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 15, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 16, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 17, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 18, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 19, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 20, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 21, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 22, ySI, 1);
       addPixel_st_mkPng(pngSTPtr, 23, ySI, 1);
   } /*Loop: pixel by pixel add*/

   addBar_st_mkPng(pngSTPtr, 0, 0, 10, 100, 3);

   print_st_mkPng(pngSTPtr, outFILE);
   errSC = 0;
   goto cleanUp_main;

   memErr_main:;
      errSC = 1;
      fprintf(stderr, "memory error%s", str_endLine);
      goto cleanUp_main;

   fileErr_main:;
      errSC = 2;
      fprintf(stderr, "file error%s", str_endLine);
      goto cleanUp_main;

   cleanUp_main:;
      if(pngSTPtr)
         freeHeap_st_mkPng(pngSTPtr);
      pngSTPtr = 0;

      if(! outFILE) ;
      else if(outFILE == stdin) ;
      else if(outFILE == stdout) ;
      else if(outFILE == stderr) ;
      else
         fclose(outFILE);
      outFILE = 0;

      return errSC;
} /*main*/
