# Use:

Here to give an idea of how to use the libraries in
  genBio. This file is/will be very long, so leverage the
  table of contents to find what you want.

The codonFun files are for converting sequences to amino
  acids or similar checks related to the codon table.

Files: codonFun.c and codonFun.h

- Dependencies:
  - ulCp.h from genLib (.h file only)
  - ntTo2Bit.h from genBio
  - revNtTo2Bit.h from genBio
  - codonTbl.h from genBio

# codonFun

Codon fun contains a set of functions to convert bases or
  amino acid sequences. It is really a wrapper around
  codonTbl, ntTo2Bit, and revNtTo2Bit.

## single amino acid conversion

The codonToAA\_codonFun function takes in three bases
  (signed characters) and returns the matching single
  letter amino acid.

```
signed char *sequenceStr = "atg";
signed char aminoAcid =
   codonToAA_codonFun(
      sequenceStr[0],
      sequenceStr[1],
      sequenceStr[2]
   );
```

The revCodonToAA_codonFun function takes in three bases
  (signed characters) and returns the matching, reverse
  complement single letter amino acid.

```
signed char *sequenceStr = "atg";
signed char reverse_complement_aminoAcid =
   revCodonToAA_codonFun(
      sequenceStr[0],
      sequenceStr[1],
      sequenceStr[2]
   );
```

## sequence translation (to amino acids)

### forward sequence translation

The seqToAA\_codonFun function translates a sequence.

| Input | Use                                     |
|:------|:----------------------------------------|
|  1st  | c-string with sequence to traslate      |
|  2nd  | c-string to hold translated sequence    |
|  3rd  | position to start translation (index 0) |
|  4th  | last base to translate                  |

Table: inputs for sequToAA\_codonFun

Make sure the 2nd input (gets ammino acid sequence) is
  large enough to hold the translated sequence
  (1 + sequence length / 3)

```
signed char *sequenceStr = "atgccctaa";
signed char *aminoAcidStr[16];
signed long lengthAASL = 0;

lengthAASL = 
   seqToAA_codonFun(
      sequenceStr,
      aminoAcidStr,
      0,
      0
   );
```

### reverse sequence translation

The revSeqToAA\_codonFun function translates a reverse
  complement sequence to an alignment.

| Input | Use                                      |
|:------|:-----------------------------------------|
|  1st  | c-string with sequence to traslate       |
|  2nd  | c-string to hold translated sequence     |
|  3rd  | first base in sequence (last translated) |
|  4th  | last base in sequence (first translated) |

Table: inputs for revSeqToAA\_codonFun

Make sure the 2nd input (gets ammino acid sequence) is
  large enough to hold the translated sequence
  (1 + sequence length / 3) or (1 + 4th input / 3). You
  can specify the from the first base in sequence (3rd
  input) to the end of the sequence by using `0` for the
  length (4th input).

```
signed char *sequenceStr = "atgccctaa";
signed char *aminoAcidStr[16];
signed long lengthAASL = 0;

lengthAASL = 
   revSeqToAA_codonFun(
      sequenceStr,
      aminoAcidStr,
      0, /*start at the frist base*/
      0  /*find the end of the sequence*/
   );
```

```
signed char *sequenceStr = "tttttttttatgccctaa";
signed char *aminoAcidStr[16];
signed long lengthAASL = 0;

lengthAASL = 
   revSeqToAA_codonFun(
      sequenceStr,
      aminoAcidStr,
      9, /*start at first a (should be 9th index)*/
      0  /*find the end of the sequence*/
   );
```

## misc functions

You can also convert a three letter amino acid to a one
  letter using aaTripToChar\_codonFun. This takes a
  c-string and returns the single letter amino acid code
  as a signed char or 0 if the input code was not an
  amino acid.

```
signed char *threeLetterStr = "Met";
signed char aaSC = aaTripToChar_codonFun(threeLetterStr);

if(! aaSC)
   /*DEAL WITH INVALID CODES*/
```

You can check if the codon is a bacterial start codon
  with codon with bactStart\_codonFun, or for reverse
  complement bactRevStart\_codonFun. The input is three
  base and the output is bacterial 1 if is a start codon
  or 0 if not.

```
signed char *seqStr = "atg";

if(bactStart_codonFun(seqStr[0], seqStr[1], seqStr[2]))
   /*have start codon*/
else
   /*not a start codon*/

if(bactRevStart_codonFun(seqStr[0], seqStr[1], seqStr[2]))
   /*have reversee complement start codon*/
else
   /*not a reverese complement start codon*/
```
