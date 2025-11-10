# Goal

Describe how to use the geneCoord struture in the
  geneCoord library. This is used to hold coordinates of
  genes.

Files: geneCoord.c geneCoord.h

- Dependenices:
  - fileFun from genLib
  - base10str from genLib
  - ulCp from genLib
  - charCp from genLib
  - genMath.h from genLib

# Overview

The geneCoord structure is designed to hold the genome
  coordinates from files from a coordinates database. It
  is by default sorted by starting position. The geneCoord
  library includes the geneCoord structure and the look
  up functions to use the structure.

# Defined variables:

There are a few pre-defined varibles that are returned or
  can affect how geneCoord works.

- def\_lenId\_geneCoord: is the number of bytes assigned
  to each gene/reference id
  - def\_lenId\_geneCoord - 1 is the maxium length allowed
    for an id/name
- def\_fileErr\_geneCoord: returned when file error
  happend
- def\_invalidEntry\_geneCoord: returned when a file had
  an invalid entry
  - might be combined with def\_fileErr\_geneCoord
- def\_memErr\_geneCoord: returned when a memory error
  happened

# geneCoord structure

The geneCoord structure has the gene name, reference name,
  start, end, and direction of genes in a tsv file. The
  length of the arrays are not stored, so you have to
  keep track of that yourselfs (weakness in my system).

- Variables:
  - idStrAry is a c-string array (maxium length 63 bytes)
    that holds the name of each gene
  - refAryStr is a c-string array (maxium length 63 bytes)
    that holds the name of each genes reference sequence
  - startAryUI has the starting positions of each gene on
    its reference
  - endAryUI has the ending positions of each gene on
    its reference
  - dirAryUC has the direction (1 reverse complement; 0
    forward) of each gene on its reference

# Flow:

The flow of the geneCoord system is read in genes from
  file, then search for genes by coordinate (start or
  range), then access values using the returned index.
  Finally, when finished free your geneCoord.

## Get and freeing gene coords

The frist step is to read in the gene coordinates file
  with the `getCoords\_geneCoord()` function. The file
  must be in the format of gene\_name, reference\_name,
  direction (`+` = foward, `-` = reverse), first
  reference coordinate, and last reference coordinate.
  The first row should be a header and is always ingored.
  The file is a tsv that can be separate by spaces and/or
  tabs.

**Example** of the input file format:

```
gene	ref	dir	start	end
atpE	NC_000962.3	+	1461045	1461290
ddn	NC_000962.3	+	3986844	3987299
eis	NC_000962.3	-	2714124	2715332
```

- Input:
  1. c-string with path to file to get coordinates from
  2. signed int pointer to get the number of genes read in
     (index 0, so is really index of the last gene)
  3. unsigned long to hold any errors
     - 0 for no errors
     - def\_fileErr\_geneCoord for file errors
     - def\_memErr\_geneCoord for memory errors
     - anything else for invalid entry error
- Output:
  - 0 for errors
  - pointer to a heap allocated geneCoord structure with
    genes in the file
    - this is sorted by starting coordinates

for the invalid entry error, you can get the line number
  of the error by doing `errorUL >> 8;` (input 3).

When finished with your geneCoord structure you can free
  it using the `freeHeap_geneCoord()` function. The input
  is the pointer to your geneCoord structure.

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "geneCoord.h"

int
main(
){
   unsigned long errorUL = 0;
   signed int lastGeneSI = 0;
   struct geneCoord *coordsHeapST =
      getCoords_geneCoord(
         (signed char *) "gene_coordinates.tsv",
         &lastGeneSI,
         &errorUL
      );
   
   if(! errorUL)
      ; /*no error*/
   else if(errorUL == def_memErr_geneCoord)
      /*deal with memory errors*/
   else if(errorUL == def_fileErr_geneCoord)
      /*deal with file errors*/
   else
   { /*Else: invalid entry in file*/
      fprintf(
         stderr,
         "line number %lu is not in the correct format\n",
         errorUL >> 8
      );
   
      /*deal with error clean up*/
   } /*Else: invalid entry in file*/
   
   if(coordsHeapST)
      freeHeap_geneCoord(coordsHeapST);
   coordsHeapST = 0;

   return 0;
} /*main*/
```

## looking up by coordinates

### looking up by gene starting coordinate

After reading in the coordinates you can then look up
  gene names by coordinates using
  the `findStart_geneCoord()` and `findRange_geneCoord()`
  functions.

The `findStart_geneCoord()` function will find a gene
  in the geneCoords structure with the input starting
  position.

- Input:
  1. geneCoord structure pointer to search for position
     (returned from getCoords\_geneCoord)
  2. starting position to search for in the geneCoord
     structure (input 1)
  3. index of last gene in the geneCoord structure
     (input 2 in getCoords\_geneCoord).
- Output:
  - index of gene that has the starting coordinate
  - -1 if the gene is not in the structure

You can then get the gene name by putting the index
  into one of the geneCoord structures arrays.

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "geneCoord.h"

int
main(
){
   unsigned long errorUL = 0;
   signed int lastGeneSI = 0;
   signed int indexSI = 0;
   struct geneCoord *coordsHeapST =
      getCoords_geneCoord(
         (signed char *) "gene_coordinates.tsv",
         &lastGeneSI,
         &errorUL
      );
   
   if(! errorUL)
      ; /*no error*/
   else if(errorUL == def_memErr_geneCoord)
      /*deal with memory errors*/
   else if(errorUL == def_fileErr_geneCoord)
      /*deal with file errors*/
   else
   { /*Else: invalid entry in file*/
      fprintf(
         stderr,
         "line number %lu is not in the correct format\n",
         errorUL >> 8
      );
   
      /*deal with error clean up*/
   } /*Else: invalid entry in file*/
   
   indexSI =
      findStart_geneCoord(coordsHeapST, 100, lastGeneSI);
   if(indexSI < 0)
      /*deal with gene is not in the array*/
   else
      fprintf(
         stdout,
         "gene %s starts at %u and ends at %u\n",
         coordsHeapST->idStrAry[indexSI],
         coordsHeapST->startAryUI[indexSI],
         coordsHeapST->endAryUI[indexSI]
      );
   
   if(coordsHeapST)
      freeHeap_geneCoord(coordsHeapST);
   coordsHeapST = 0;

   return 0;
} /*main*/
```

### looking up by range

You can also look up the gene(s) that may match a range
  of coordinates by using the `findRange_geneCoord()`
  function. this will return the index of the frist gene
  that can fall in your range.

- Input:
  1. geneCoord structure pointer to search for position
     (returned from getCoords\_geneCoord)
  2. starting position to search for in the geneCoord
     structure (input 1)
  3. ending position to search for in the geneCoord
     structure (input 1)
  4. index of last gene in the geneCoord structure
     (input 2 in getCoords\_geneCoord).
- Output:
  - first gene that is in the range
  - -1 if the gene is not in the structure

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "geneCoord.h"

int
main(
){
   unsigned long errorUL = 0;
   signed int lastGeneSI = 0;
   signed int indexSI = 0;
   struct geneCoord *coordsHeapST =
      getCoords_geneCoord(
         (signed char *) "gene_coordinates.tsv",
         &lastGeneSI,
         &errorUL
      );
   
   if(! errorUL)
      ; /*no error*/
   else if(errorUL == def_memErr_geneCoord)
      /*deal with memory errors*/
   else if(errorUL == def_fileErr_geneCoord)
      /*deal with file errors*/
   else
   { /*Else: invalid entry in file*/
      fprintf(
         stderr,
         "line number %lu is not in the correct format\n",
         errorUL >> 8
      );
   
      /*deal with error clean up*/
   } /*Else: invalid entry in file*/
   
   indexSI =
     findRange_geneCoord(coordsHeapST, 100, 200, lastGeneSI);
   if(indexSI < 0)
      ; /*no geens were in this range*/
   else
   { /*Else: found genes in range*/
      while(coordsHeapST->startAryUI[indexSI] <= 200)
      { /*Loop: print genes that were in range*/
         fprintf(
            stdout,
            "gene %s starts at %u and ends at %u\n",
            coordsHeapST->idStrAry[indexSI],
            coordsHeapST->startAryUI[indexSI],
            coordsHeapST->endAryUI[indexSI]
         );
   
         ++indexSI;
      } /*Loop: print genes that were in range*/
   } /*Else: found genes in range*/
   
   if(coordsHeapST)
      freeHeap_geneCoord(coordsHeapST);
   coordsHeapST = 0;

   return 0;
} /*main*/
```

## looking up by name

You can look up genes by name, but you first need to sort
  the arrays in the geneCoord structure by name using
  the `sortName_geneCoord();` function. You can resort
  by starting position using the `sort_geneCoord()`
  function.

After sorting you can use the `findName_geneCoord()`
  function to look up genes by name.

### sorting by name

You can sort the geneCoordinates structure by name using
  the `sortName_geneCoord()` function.

- Input:
  1. geneCoord structure pointer to search for position
     (returned from getCoords\_geneCoord)
  2. 0 to specify the first gene
  3. index of last gene in the geneCoord structure
     (input 1)

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "geneCoord.h"

int
main(
){
   unsigned long errorUL = 0;
   signed int lastGeneSI = 0;
   signed int indexSI = 0;
   struct geneCoord *coordsHeapST =
      getCoords_geneCoord(
         (signed char *) "gene_coordinates.tsv",
         &lastGeneSI,
         &errorUL
      );
   
   if(! errorUL)
      ; /*no error*/
   else if(errorUL == def_memErr_geneCoord)
      /*deal with memory errors*/
   else if(errorUL == def_fileErr_geneCoord)
      /*deal with file errors*/
   else
   { /*Else: invalid entry in file*/
      fprintf(
         stderr,
         "line number %lu is not in the correct format\n",
         errorUL >> 8
      );
   
      /*deal with error clean up*/
   } /*Else: invalid entry in file*/
   
   sortName_geneCoord(coordsHeapST, 0,lastGeneSI);
   
   if(coordsHeapST)
      freeHeap_geneCoord(coordsHeapST);
   coordsHeapST = 0;

   return 0;
} /*main*/
```

### sorting by starting position

In some cases you may want to undo the sorting you did by
  gene name. In those cases use the `sort_geneCoord()`
  function.

- Input:
  1. geneCoord structure pointer to search for position
     (returned from getCoords\_geneCoord)
  2. 0 to specify the first gene
  3. index of last gene in the geneCoord structure
     (input 1)

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "geneCoord.h"

int
main(
){
   unsigned long errorUL = 0;
   signed int lastGeneSI = 0;
   signed int indexSI = 0;
   struct geneCoord *coordsHeapST =
      getCoords_geneCoord(
         (signed char *) "gene_coordinates.tsv",
         &lastGeneSI,
         &errorUL
      );
   
   if(! errorUL)
      ; /*no error*/
   else if(errorUL == def_memErr_geneCoord)
      /*deal with memory errors*/
   else if(errorUL == def_fileErr_geneCoord)
      /*deal with file errors*/
   else
   { /*Else: invalid entry in file*/
      fprintf(
         stderr,
         "line number %lu is not in the correct format\n",
         errorUL >> 8
      );
   
      /*deal with error clean up*/
   } /*Else: invalid entry in file*/
   
   sortName_geneCoord(coordsHeapST, 0, lastGeneSI);
   sort_geneCoord(coordsHeapST, 0, lastGeneSI);
   
   if(coordsHeapST)
      freeHeap_geneCoord(coordsHeapST);
   coordsHeapST = 0;

   return 0;
} /*main*/
```

### Looking up gene by name

After sorting by gene name you can search for gene names
  in a geneCoord structure using
  the `findName_geneCoord()` function.

- Input:
  1. geneCoord structure pointer to search for position
     (returned from getCoords\_geneCoord)
  2. c-string with name of gene to search for
  3. index of last gene in the geneCoord structure
     (input 1)
- Output:
  - index of matching gene (if found)
  - -1 if gene name was not in the coordinates structure

**Example**

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "geneCoord.h"

int
main(
){
   unsigned long errorUL = 0;
   signed int lastGeneSI = 0;
   signed int indexSI = 0;
   struct geneCoord *coordsHeapST =
      getCoords_geneCoord(
         (signed char *) "gene_coordinates.tsv",
         &lastGeneSI,
         &errorUL
      );
   
   if(! errorUL)
      ; /*no error*/
   else if(errorUL == def_memErr_geneCoord)
      /*deal with memory errors*/
   else if(errorUL == def_fileErr_geneCoord)
      /*deal with file errors*/
   else
   { /*Else: invalid entry in file*/
      fprintf(
         stderr,
         "line number %lu is not in the correct format\n",
         errorUL >> 8
      );
   
      /*deal with error clean up*/
   } /*Else: invalid entry in file*/
   
   sortName_geneCoord(coordsHeapST, 0, lastGeneSI);
   
   indexSI =
      findName_geneCoord(
         coordsHeapST,
         (signed char *) "eis",
         lastGeneSI
      );
   
   if(indexSI < 0)
      ; /*eis gene is not in the coordinates file*/
   else
      fprintf(
         stdout,
         "eis starts at %u and ends at %u\n",
         coordsHeapST->idStrAry[indexSI],
         coordsHeapST->startAryUI[indexSI],
         coordsHeapST->endAryUI[indexSI]
      );
   
   if(coordsHeapST)
      freeHeap_geneCoord(coordsHeapST);
   coordsHeapST = 0;

   return 0;
} /*main*/
```
