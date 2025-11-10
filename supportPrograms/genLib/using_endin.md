# Goal:

Descripe how to convert endin of a number using the endin
  library.

Files: endin.h and endin.c

# endin

Converts integers (or longs/shorts) to big endin or little
  endin format. Not very great, but should allow
  variations in datatype size.

For big endin OS's use `-DBIG_ENDIN` flag.

- big endin Conversion functions
  - these are only enabled if you do not use `-DBIG_ENDIN`
    otherwise the input is returned
  - long: bigEndUL = ulToBig\_endin(input)
  - int: bigEndUI = uiToBig\_endin(input)
  - short: bigEndUS = usToBig\_endin(input)
- little endin Conversion functions
  - these are only enabled if `-DBIG_ENDIN` is used,
    otherwise the input is returned
  - long: littleEndUL = ulToLittle\_endin(input)
  - int: littleEndUI = uiToLittle\_endin(input)
  - short: littleEndUS = usToLittle\_endin(input)
- swap endin
  - these always convert to the oppisiste endin
  - long: swapEndUL = ulSwap\_endin(input)
  - int: swapEndUI = uiSwap\_endin(input)
  - short: swapEndUS = usSwap\_endin(input)

**Example**

```
unsigned long valueUL = 10;
unsigned long bigEndinValueUL = ulToBig_endin(valueUL);
   /*this will make sure valueUL is in bigendin format*/
```
