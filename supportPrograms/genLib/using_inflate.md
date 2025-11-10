# Goal:

Describe how to use the inflate library to uncompress
  gz files.

Files: inflate.h and inflate.c

- Dependencies:
  - genMath
  - ulCp
  - endin
  - checkSum
  - endLine.h

# inflate

I would recomend using GNU's .gz program or zlib instead.
  Only use this if they do not provide what you need (
  very unlikely).

Maybe not best setup, but should work on most OS's. For
  16 bit (not recomemded and untested) OS's you will need
  at least 40kb to 50kb of memmory.

For inflate you will need to work with the file\_inflate
  (st02 in inflate.h) structure. This structure manages
  the buffers needed to read in the zip file and also
  store and uncompress the data.

Most of the variables should not be interacted with,
  however, there are a few that are safe.

- timeAryUC is a four byte unsinged char array with the 
  date stamp from a gz header (currently one). You can
  cast this to an int
- nameStr is a c-string with the file name
  - nameLenUC has the name length
- buffStr is the uncompressed buffer. Do not fiddle with
  this, unless you know what you are doing
  - buffPosSL is the current position in the uncompressed
    buffer

| Flag            | Use                             |
|:----------------|:--------------------------------|
| `-DNOUL`        | see ulCp                        |
| `-DSIXTEEN_BIT` | do not use gz (32 bit) checksum |

Table: flags for inflate. For Windows replace `-D`
  with `/D`

## file\_inflate overview

To use a file\_inflate struct, you must first initialize
  it, then add your file to it, then set the file you
  added to null (inflate will call fclose), then read
  bytes from the buffer. When finished you should free
  your structure.

## file\_inflate initialize

To initialize use init\_file_\inflate (fun08 inflate.c/h).
  This takes in a pointer to an inflate structure to
  initialize.

## file\_inflate add file

To add a file call add\_file\_inflate (fun25 inflate.c/h).
  This takes in a FILE pionter (file to add), a
  pointer to a file\_inflate structure to add the file to,
  and a null (or two byte char array with the first two
  bytes of file). It returns 0 for success and a wide
  variarity of errors (see function header for a list)
  for problems. After calling, you should set the FILE
  pointer you used to 0 (inflate will close the file with
  fclose).

## file\_inflate get data

You can use get\_inflate (fun23 inflate.c/h) to read in
  uncompressed bytes by count (like fread), by line
  (like fget), or by block with get\_inflate. The input
  is:

1. a pointer to a file\_inflate structure to inflate
2. the maximum number of bytes to read in
  - depnding on choosen method it could be less
3. a flag (endBl) marking if you want to behave like
  fread (0) (maximum bytes or end of deflate), fget
  (1) (line by line), or read block by block (2)
  - this is limited to the maximum number of bytes to
    grab, so it may not always get the entire line or
    block
  - use `4` to tell inflate there are always two extra
    btyes in the buffer. This sets inflate to grab line
    breaks when buffer seems full.
  - use `1 | 4` to apply `4` when just grabbing one line
  - use `2 | 4` to apply `4` when grabbing a block
4. the c-string buffer to save the uncompressed data to
  - on text based compression, line endings are changed to
    OS line endings
5. a pointer to a signed char to hold any errors
  - set to def\_eof\_inflate if finished the final block
  - see the codumentation

The output is the number bytes read in or a negative
  number for errors. Also the input c-string buffer (4)
  will be filled with the uncompressed data.

## file\_inflate checksum and move on

After finishing an deflate entry, you can move to the next
  deflate entry (and do checksum) with next\_file\_inflate
  (fun24 inflate.c/h). The input is:

- a ponter to a file\_inflate structure to move to the next
  entry on
- 0; here for add\_file\_inflate
  - would be next to bytes (if already read in)
- 0; here for add\_file\_inflate
  - 1: marks the first time reading this file

The return value is 0 for no errors, def\_eof\_inflate if
  at the end of the file, and several other return values
  for errors. Check the function header for a list.

## file\_inflate next file or re-read

You can reset the inflate structure by using
  blank\_file\_inflate. This will either close the file or
  seek the start of the file. It will also set all
  variables to 0.

- Input:
  - pointer to file\_inflate structure to blank
  - 1 to close (fclose) the file and set to 0
    - 0 to go to the start of the file

This funcition has no return value.

## file\_inflate clean up

You can clean up your file\_inflate structure using
  freeStack\_file\_inflate to free the variables in a
  file\_inflate (fun09 inflate.c/h) structure or
  freeHeap\_file\_inflate (fun10 inflate.c/h) to free
  a file\_inflate structure (you must set pointer to
  0/null).

## Inflate example

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>
#include "../genLib/inflate.h"

/*inflate also uses ulCp.h endin.h, checksum.h, genMath.h,
`   and endline.h from genLib
*/

int main(){
   signed char errSC = 0;
   signed char outStr[1024 + 1];
   signed long bytesSL = 0;

   struct file_inflate zipST;
   FILE *zipFILE = stdin;

   init_file_inflate(&zipST);

   if(add_file_inflate(zipFILE, &zipST, 0))
      /*HANDEL ERRORS*/
   zipFILE = 0; /*so I do not use fclose at end*/
   
   while(! errSC)
   { /*Loop: ucompress and print file*/
      uncompress:;

      bytesSL =
         get_inflate(
            &zipST,
            1024 - 2,/*-2 to add 2 extra bytes at end*/
            4,      /*grab line break if have 1022 chars*/
            outStr,  /*buffer to hold new bytes*/
            &errSC   /*error messages*/
         );
         /*change 4 to `1 | 4` = 5 to read line by line*/

      if(bytesSL < 0)
         /*handel errors*/

      printf("%s", outStr);
   } /*Loop: ucompress and print file*/

   errSC = next_file_inflate(&zipST, 0, 0);

   if(! errSC)
      goto uncompress; /*restart uncompress loop*/
   else if(errSC != def_eof_inflate)
      /*HANDEL ERRORS*/

   freeStack_file_inflate(&zipST);
      /*will use fclose on the FILE pointer in zipST*/
   return 1;
}
```
