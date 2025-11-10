# Goal:

Describe purpose and use of the indexToCoord library.

Files: indexToCoord.c and indexToCoord.h

# Overview:

The index to coord library is designed to conver the
  indexSL variable in a dirMatrix structure into a
  reference, query, or reference and query coordinates.
  This is mainly for Needlmans and Watermans, so you
  should not have to mess with this library.

There are three functions, qryCoord\_indexToCoord,
  refCoord\_indexToCoord, and indexToCoord.

For the qryCoord\_indexToCoord and refCoord\_indexToCoord
  functions take the same input. The qry returns the
  query coordinate, wile the ref returns the reference
  coordinate.

- Input (qry and ref):
  1. length of reference (refLenSL in the dirMatrix
     structure)
  2. coordinate to convert (indexSL in the dirMatrix
     structure).

```
unsigned long queryUL =
   qryCoord_indexToCoord(<reference_length>,<coordinate>);
unsigned long referenceUL =
   refCoord_indexToCoord(<reference_length>,<coordinate>);
```

The indexToCoord function modifies two input long pointers
  to hold the query and the reference coordinate. Think of
  it as refCoord\_indexToCoord and qryCoord\_indexToCoord
  combined together.

- Input:
  1. length of reference (refLenSL in the dirMatrix
     structure)
  2. coordinate to convert (indexSL in the dirMatrix
     structure).
  3. unsigned long pointer to get the reference coordiante
  4. unsigned long pointer to get the query coordiante

```
unsigned long queryUL = 0;
unsigned long referenceUL = 0;

indexToCoord(
   <reference_length>,
   <coordinate>,
   &referenceUL, /*gets reference coordinate*/
   &queryUL,     /*gets query coordinate*/
);
```
