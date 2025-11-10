# Use:

Here to give you an orientation of the C code files and
  functions. This is here to list what is present and
  what each library does. To see more specific details
  check the `using_<library>.md` files.

# Libraries:

- Table Of Contents (TOC)
  - base10str (see using\_base10str.md)
    - c-string to base10 number conversoin
  - charCp (see using\_charCp.md)
    - byte copy functions
  - checkSum (see using\_checkSum.md)
    - does checks sums for inflate; currently crc32 and
      alder32 only
  - endin (see using\_endin.md)
    - converts datatypes to bigendin format
  - endLine (see using\_endLine.md)
    - line ending variable for OS's
  - fileFun (see using\_fileFun.md)
    - randon file funcitions, currently only a count
      number lines functions
  - genMath (see using\_genMath.md)
    - math (branchless) functions
  - inflate (see using\_inflate.md)
    - inflates deflated data (for .gz files)
  - mkPng (see using\_mkPng.md)
    - not working; goal is make uncompressed pallete pngs
    - documentation is not complete, will finish when I
      get mkPng working
  - numToStr (see using\_numToStr.md)
    - number to c-string
  - ptrAry (see using\_ptrAry)
    - c-string array (using pointers)
  - shellSort (see using\_shellShort)
    - shell short and their functions
  - strAry (see using\_strAry)
    - c-string array that is in a (signed char \*) array,
      so no structure, but maximum string lengths are
      63 and will waste space on short strings.
  - ulCp (see using\_ulCp)
    - copying strings using unsigned longs (8 bytes at a
      time). This is slower then system functions.
