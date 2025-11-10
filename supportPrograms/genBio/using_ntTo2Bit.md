# Use:

Here to give an idea of how to use the ntTo2Bit look up
  array in genBio.

Look up table to convert nucleotides to two bit numbers.
  A sepecial 3rd bit is used to store errors (non ATCG).


Files: ntTo2Bit.h

Dependencies: None

# ntTo2Bit

Is a lookup table to convert nucleotides to array lookup
  index's. The lookup table is called ntTo2Bit. It takes
  in a unsigned character and returns the converted
  nucleotide back as an unsigned character.

| Nucleotide  | Return                   |
|:------------|:-------------------------|
| T           | def\_t\_ntTo2Bit         |
| C           | def\_c\_ntTo2Bit         |
| A           | def\_a\_ntTo2Bit         |
| G           | def\_g\_ntTo2Bit         |
| other input | def\_err3rdBit\_ntTo2Bit |

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
