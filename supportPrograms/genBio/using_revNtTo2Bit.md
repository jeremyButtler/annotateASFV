# Use:

Here to give an idea of how to use the revNtTo2Bit look
  up table in genBio.

Look up table to first complement and then convert
  nucleotides to two bit numbers.  A sepecial 3rd bit is
  used to store errors (non ATCG).

Files: revNtTo2Bit.h

Dependencies: None

# revNtTo2Bit

Is a lookup table to reverese complement and convert
  nucleotides to array lookup index's. The lookup table is
  called revNtTo2Bit. It takes in a unsigned character and
  returns the converted nucleotide back as an unsigned
  character.

| Nucleotide  | Return                      |
|:------------|:----------------------------|
| T           | def\_a\_revNtTo2Bit         |
| C           | def\_g\_revNtTo2Bit         |
| A           | def\_t\_revNtTo2Bit         |
| G           | def\_c\_revNtTo2Bit         |
| other input | def\_err3rdBit\_revNtTo2Bit |

Simple example; printing a translated, reverse complement
  sequence:

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
   unsigned char indexUC = 5;
   unsigned char ntAryUC[16];
   unsigned char indexUC = 0;
   
   while(sequenceStr[ntUC])
      ntAryUC[indexUC++] =
         revNtTo2Bit[(unsigned char) sequenceStr[ntUC++]];
         /*filling in backwards to reverese sequence*/

   indexAryUC[ntUC] = '\0';
   indexAryUC[ntUC + 1] = '\0';
   indexAryUC[ntUC + 2] = '\0';

   ntUC = 0;

   while(indexAryUC[ntUC])
   { /*Loop: convert index to amino acids*/
      putc(
         stdout,
         codonTbl[
              ntAryUC[ntUC]     /*part of codonTbl line*/
            ][ntAryUC[ntUC + 1] /*part of codonTbl line*/
            ][ntAryUC[ntUC + 2] /*part of codonTbl line*/
         ] /*end of codonTbl line*/
      );
      ntUC += 3; /*move to next amino acid
   } /*Loop: convert index to amino acids*/
}
```
