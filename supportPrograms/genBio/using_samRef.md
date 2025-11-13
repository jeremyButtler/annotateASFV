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
  find/extract reference ids and lengths. This also can
  be used to bin reads in a sam file by mapped reference.

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
  
## refs\_samRef workflow

### overview

The workflow for a refs\_samRef struct is similar to a
  samEntry struct. You first initialize the struct, then
  setup initial memory. After that you can read in a
  samEntry header or add references one by one. You can
  then search the structure for references. Finally you
  will free the struct.

### setup

You can initialize a refs\_samRef struct with
  `init_refs_samRef(refs_samRefPointer);`.

- You can then setup memory using
  - Note: as of 2025-11-12, the getRefLen\_samRef and
    addRef\_samRef functions will check if things need
    to be setup

```
init_refs_samRef(refs_samRefPointer);
if(setup_refs_samRef(refs_samRefPointer);
   /*deal with memory errors*/
```

### add references

You can then read the headers from a sam file using
  getRefLen\_samRef. This will read the headers until
  it hits the frist sequence.

- Input:
  1. pointer to a `refs\_samRef` structure to store
     lengths and reference ids in
  2. pointer to a samEntry to get lines
     - will be set to first sequence line
  3. FILE pointer to sam file to get references from
  4. FILE pointer to file to print headers to
     - use 0/null to not print headers
  5. c-string pointer to get non-reference header entries
     - resized as needed (ok if pointer is set to 0)
  6. unsinged long pointer to get/has header c-string size
- Output:
  - Returns 0 for no errors
  - Returns def\_memErr\_samEntry for memory errors
  - Returns def\_fileErr\_samEntry for file errors
  - stores header in input c-string
  - updates unsigned long pointer to header c-string is
    resized

Another option is to add each reference inividually with
  the addRef_samEntry function. This will add in the
  reference id to the correct (sorted) position.

- Input:
  1. c-string with reference id/name to add
  2. length of the reference genome
  3. refs_samRef structer pointer to add the reference to
  4. signed char pointer to get error values
- Output
  - Returns: -1 for errors
  - Returns: index of the inserted reference
    - note, this could change by one position for each new
      reference you add, so it is temporary
  - Modifies: errSCPtr to be 0 for no errors
  - Modifies: errSCPtr to be def_expand_samEtnry if needed
    to resize the memory usage
    - you can ignore this error
  - Modifies; errSCPtr to be def_memErr_samEntry if had
    memory errors

### finding ids

You can then find references in the refs\_samRef struct
  using findRef\_refs\_samRef. The input is the id to
  search for and teh refs\_samRef struct to search. The
  output is the index of the reference or a negative
  number if no referene was found.

### binning reads

For read binning you have two choices, you can bin reads
  into one bin for each reference sequence or you can
  have samRef find potential genbank assembly accesions
  and try to merge the accesions.

**find bins**

For merging assembly reference sequences you first need
  to create the bins. You can do this with
  buildRefMergeIndex\_samRef.

- Input:
  1. refs_samRef structure pionter with reference ids to
     search for merging
- Output:
  - signed int array with the index to assign the id to
    - the length of the array is the number of ids in
      input 1
    - values:
      - -1 if reference id belongs to its own bin
      - index of reference id's bin this id belongs to

**print reads to bins**

You can then print the reads to bins using the
  pSamToRef\_samRef function. This will open a file for
  you and print the read and if needed header to the
  output file. This is a little slow because of all the
  file opening, but it does work.

- Input:
  1. samEntry structure pionter with read to print
  2. c-string with the prefix for the output file name
     - `<prefix>-<reference_id>.sam`
  3. c-string with the general header (no reference
     information) to print
     - input 5 of getRefLen\_samRef
     - `@HD` entry is printed first
     - all other entries are printed after reference ids
  4. 1 to print unmapped reads to the unmapped file or
     0 to not print unmapped reads 
  5. refs_samRef structure pointer with reference ids
  6. signed int array or 0/null with the indexs for each
     references bin
     - return value from buildRefMergeIndex_samRef
     - if 0/null, each reference sequence gets its own bin
- Output
  - Prints: read and if needed header to the bin the read
    was assigned to
  - Returns: 0 for success
  - Returns: 1 if read did not have a reference id
  - Returns: 2 if read had no sequence
  - Returns: -1 for file errors

**printing headers (not reccomeneded)**

I would not recomend using this function, because the
  pSamToRef\_samRef fuction does this for you, but you can
  print the header for each bin sam file using the
  pSamHeader_samRef function. The only reason to use this
  would be if you were building your own binning function.

- Input:
  1. index of the reference to print header for
  2. general header (no reference information) to print
     - `@HD` entry is printed first
     - all other entries are printed after reference ids
  3. 1 to print unmapped reads to the unmapped file or
     0 to not print unmapped reads 
  4. refs_samRef structure pointer with reference ids
  5. signed int array or 0/null with index of bins to
     print out
     - return value from buildRefMergeIndex_samRef
  6. FILE pointer to print the header to
- Output
  - prints header to input 6

### clean up

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

Also, if you called buildRefMergeIndex\_samRef, remember
  to free the returned signed int array with `free()`.

## refs\_samRef example

### getting reference lengths

```
#ifdef PLAN9
  #include <u.h>
  #include <libc.h>
#else
  #include <stdlib.h>
#endif

#include <stdio.h>
#include "samEntry.h"
#include "samRef.h"


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

### binning reads

```
#ifdef PLAN9
  #include <u.h>
  #include <libc.h>
#else
  #include <stdlib.h>
#endif

#include <stdio.h>
#include "samEntry.h"
#include "samRef.h"


int
main(
){
   signed char errorSC = 0;
   signed char *headerStr = 0;
   unsigned long headerSizeUL = 0;
   signed int *indexHeapArySI = 0;
   
   struct refs_samRef refsStruct;
   struct samEntry samStruct;
   
   FILE *samFILE = fopen("file.sam", "r");
   
   init_refs_samRef(&refsStruct);
   init_samEntry(&samStruct);
   
   if(! samFILE)
      /*deal with file errors*/
   
   errorSC =
      getRefLen_samRef(
         &refStruct,
         &samStruct, /*set to first read in sam file*/
         samFILE,
         0,      /*not printing headers to output file*/
         &headerStr,
         &headerSizeUL
      );
   
   indexHeapArySI =
      buildRefMergeIndex_samRef(&refsStruct);
      /*you can skip this step to have each reference
      `  id have its own bin
      */

   while(! errorSC)
   { /*Loop: bin read*/
      if(samStruct.flagUS & (4 | 256 | 2048))
         goto nextRead_main;

      errorSC =
         pSamToRef_samRef(
            &samStruct,            /*has read to print*/
            (signed char *) "out", /*name of file*/
            headerStr,             /*header of file*/
            0,                 /*discared unmapped reads*/
            &refStruct,            /*has reference ids*/
            indexHeapArySI         /*has bins to merge*/
         );

      if(errorSC < 0)
         /*deal with file errors*/

      nextRead_main:;
         errorSC = get_samEntry(&samStruct, samFILE);
   } /*Loop: bin read*/

   if(errorSC != def_EOF_samEntry)
      /*deal with error reading the sam file*/

   freeStack_samEntry(&samStruct);
   freeStack_refs_samRef(&refStruct);
   fclose(samFILE);
   samFILE = 0;
   
   if(headerStr)
      free(headerStr);
   headerStr = 0;

   if(indexHeapArySI)
      free(indexHeapArySI);
   indexHeapArySI = 0;

   return 0;
} /*main*/
```
