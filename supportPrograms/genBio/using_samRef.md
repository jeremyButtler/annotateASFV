# Use:

Here to give an idea of how to use the samRef library in
  genBio.

Functions for extracing reference information from a
  sam file. Uses the samEntry system.

Files: samRef.c and samRef.h

- Dependencies:
  - samEntry from genBio
  - strAry from genLib
  - ulCp from genLib
  - base10str from genLib
  - endLine.h from genLib

- Hidden dependencies:
  - numToStr from genLib
  - fileFun from genLib
  - ntTo5Bit from genBio

You should be familar with the samEntry struct (see
  using\_samEntry.md) and strAry (see
  ../genLib/using\_genLib.h) before using samRef.

# samRef

Has structure and function to scan sam file header and
  find/extract reference ids and lengths.

# refs\_samRef:

Has the list of references and their lengths found in the
  sam file header.

- lenAryUI array that has all the reference lengths
- idAryStr has the reference names
  - this is a c-string array from ../genLib/strAry.c, so
    you will need to use get\_strAry() (fun03
    ../genLib/strAry.c) to get ids out of it.
- numRefUI has the number of references in lenAryUI and
  idAryStr
  
### refs\_samRef workflow

The workflow for a refs\_samRef struct is similar to a
  samEntry struct. You first initialize the struct, then
  setup initial memory. After that you can read in a
  samEntry header or add references one by one. You can
  then search the structure for references. Finally you
  will free the struct.

You can initialize a refs\_samRef struct with
  `init_refs_samRef(refs_samRefPointer);`.

You can then setup memory using 

```
init_refs_samRef(refs_samRefPointer);
if(setup_refs_samRef(refs_samRefPointer);
   /*deal with memory errors*/
```

You can then read the headers from a sam file using
  getRefLen\_samEntry. This will read the headers until
  it hits the frist sequence.

- Input:
  - pointer to a `refs\_samRef` structure to store
    lengths and reference ids in
  - pointer to a samEntry to get lines
    - will be set to first sequence line
  - FILE pointer to sam file to get references from
  - FILE pointer to file to print headers to
    - use 0/null to not print headers
  - c-string pointer to get non-reference header entries
    - resized as needed (ok if pointer is set to 0)
  - unsinged long pointer to get/has header c-string size

- Output:
  - Returns 0 for no errors
  - Returns def\_memErr\_samEntry for memory errors
  - Returns def\_fileErr\_samEntry for file errors
  - stores header in input c-string
  - updates unsigned long pointer to header c-string is
    resized

You can then find references in the refs\_samRef struct
  using findRef\_refs\_samRef. The input is the id to
  search for and teh refs\_samRef struct to search. The
  output is the index of the reference or a negative
  number if no referene was found.

Finally you can free the refs\_samRef structure using
  freeStack\_refs\_samRef (for variables in struct) or
  with freeHeap\_refs\_samRef (for entire struct). If
  you use freeHeap\_refs\_samRef, remember to set you
  structure pointer to null.

```
freeStack_refs_samRef(refs_samRefPointer);
```

```
freeHeap_refs_samRef(refs_samRefPointer);
refs_samRefPointer = 0;
```

### refs\_samRef example

```
#ifdef PLAN9
  #include <u.h>
  #include <libc.h>
#else
  #include <stdlib.h>
#endif

#include <stdio.h>
#include "samEntry.h"


int
main(
){

   signed char errorSC = 0;
   signed char *headerStr = 0;
   unsigned long headerSizeUL = 0;
   signed long indexSL = 0;
   
   struct refs_samRef refsStruct;
   struct samEntry samStruct;
   
   FILE *samFILE = fopen("file.sam", "r");
   
   init_refs_samRef(&refsStruct);
   init_samEntry(&samStruct);
   
   if(! samFILE)
      /*deal with file errors*/
   if(setup_samEntry(&samStruct)
      /*deal with memory errors*/
   if(setup_refs_samRef(&refsStruct))
      /*deal with memory errors*/
   
   errorSC =
      getRefLen_samRef(
         &refStruct,
         &samStruct, /*set to first read in sam file*/
         samFILE,
         0,      /*not printing headers to output file*/
         &headerStr,
         &headerSizeUL
      );
   
   indexSL =
     findRef_refs_samRef(&refStruct,samStruct.refIdStr);
   
   if(indexSL < 0)
      /*likely umapped read (reference not in header)*/
   else
   { /*Else: found reference, print reference length*/
      printf(
         "%s\t%li\n",
         refStruct.idAryStr[indexSL],
         refStruct.lenAryUI[indexSL]
      );
   } /*Else: found reference, print reference length*/
   
   freeStack_samEntry(&samStruct);
   freeStack_refs_samRef(&refStruct);
   fclose(samFILE);
   samFILE = 0;
   
   if(headerStr)
      free(headerStr);
   headerStr = 0;
   
   return 0;
}
```
