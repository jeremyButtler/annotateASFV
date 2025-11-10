# Use:

Here to give an idea of how to use the codonTbl header in
  genBio.

The codonTbl file is the look up array to convert
  nucleotides to the correct codon. Before using you must
  the nucleotides using the ntTo2Bit or revNtTo2Bit look
  up table.

Files: codonTbl.h

Dependencies: None, but ntTo2Bit or revNtTo2Bit is
  recommended.

# codonTbl

The codon table is a lookup table to convert three bases
  into an amino acid. It uses single letter amino acid
  symbols. There is a list of the symbols and a codon
  table at the bottom of the file. The tables name is
  codonTbl and it takes in three index (bases). You can
  convert an ascii base to a lookup value using ntTo2Bit
  or revNtTo2Bit. The return value is the ammino acid
  letter (x for any) as a signed char.

See using\_ntTo2Bit.md or using\_revNtTo2Bit.md for an
  example of how to use codonTbl.

Simple example; printing a translated sequence:

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include <stdio.h>
#include "ntTo2Bit.h"
#include "codonTbl.h"

int
main(
){
   signed char *sequenceStr = "atgccc";
   unsigned char ntUC = 0;
   unsigned char ntAryUC[16];
   unsigned char indexUC = 0;
   
   while(sequenceStr[ntUC])
      ntAryUC[ntUC] =
         ntTo2Bit[ (unsigned char) sequenceStr[ntUC++] ];

   indexAryUC[ntUC] = '\0';
   indexAryUC[ntUC + 1] = '\0';
   indexAryUC[ntUC + 2] = '\0';

   ntUC = 0;

   while(indexAryUC[ntUC])
   { /*Loop: convert index to amino acids*/
      putc(
         stdout,
         codonTbl[
              ntAryUC[ntUC]      /*part of codonTbl line*/
            ][ntAryUC[ntUC + 1]  /*part of codonTbl line*/
            ][ntAryUC[ntUC + 2]  /*part of codonTbl line*/
         ] /*end of codonTbl line*/
      );
      ntUC += 3; /*move to next amino acid
   } /*Loop: convert index to amino acids*/
}
```
