# Use:

Here to give an idea of how to use the ntTo5Bit look up
  array in genBio. Not the best system.


Files: ntTo5Bit.h

Dependencies: None

# ntTo5Bit

Look up table to convert nucleotides to a five bit number.
  The first four bits record the nucleotide, while the
  fifth bit records if it was an anonymous nucleotide.
  An extra 6th bit is used to tell if there was an
  error (1) or no errors (0).

This is not a very efficent method for space.

| Nucleotide  | Return                         |
|:------------|:-------------------------------|
| T           | def\_t\_four\_ntTo2Bit         |
| C           | def\_c\_four\_ntTo2Bit         |
| A           | def\_a\_four\_ntTo2Bit         |
| G           | def\_g\_four\_ntTo2Bit         |
| anonymous   | def\_n\_fithBit\_ntTo5Bit      |
| error       | def\_err\_sixBit\_ntTo5Bit     |

Simple example; printing a translated sequence:

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#endif

#include "ntTo5Bit.h"

int
main(
){
   signed char *sequenceStr = "atgccc";
   unsigned char ntUC = 0;
   unsigned char ntAryUC[16];
   unsigned char indexUC = 0;
   
   while(sequenceStr[ntUC])
   { /*Loop: get nucleotide codes*/
      ntAryUC[ntUC] =
         ntTo5Bit[ (unsigned char) sequenceStr[ntUC++] ];

      if(ntAryUC[ntUC] & def_err_sixBit_ntTo5Bit)
         /*not a valid base*/
      else if(ntAryUC[ntUC] & def_n_fifthBit_ntTo5Bit)
         /*is an anonymous base*/
      else
         /*else is a regular base*/
   } /*Loop: get nucleotide codes*/

   /*code to use converted sequence goes here*/


  return 0;
} /*main*/
```
