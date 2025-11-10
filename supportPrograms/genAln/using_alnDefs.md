# Goal:

Give a brief overview on how to use the alnDefs defaults
  to change settings for the alingment programs

Files: alnDefs.h

# alnDefs

This has the default settings for the Needleman, Waterman,
  and Hirschberg aligners and the functions that use them.
  The main thing of note is that it as the gap scoring
  penalties and the scoring matrix.

Do not mess around with the `def_mv*_alnDefs` defined
  variables or any other variable in Sec02. These are used
  for the directional matrix.

- In Sec03: you have the gap scoring variables.
  - def\_gapOpen\_alnDefs is the gap opening penalty
    - should be negative
  - def\_gapExtend\_alnDefs is the gap extension penalty
    - should be negative
  - def\_scoreAdj\_alnDefs is what to dived the final
    score by before reporting

- Sec04 has the scoring matrix. All it is, is a list of
  defined variabls. One for each nucleotide pairing.
  - Format: `def_<query>To<reference>_alnDefs`
    - `<query>` is the first nucleotide to compare
    - `<reference>` is the second nucleotide to compare
    - Example: for A compared to T is `def_AToT_alnDefs`
  - The number assigned is the score for that pair

- Sec05 has the match matrix
  - Format: `def_<query>Eql<reference>_alnDefs`
    - `<query>` is the first nucleotide to compare
    - `<reference>` is the second nucleotide to compare
  - Possible values:
    - 0 for no match
    - def\_ntEql\_alnDefs for a match
    - def\_anonymous\_alnDefs if match, but it had an
      anonymous base (W,S,M,K,R,Y,B,D,H,V,N,X)
      - X is for amino acids

