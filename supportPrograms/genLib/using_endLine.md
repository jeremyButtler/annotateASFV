# Goal:

Describe how to use the `str_endLine` variable in the
  `endLine` library to print line endings.

Files: endLine.h

# endLine

Has the line endings for each OS. Only has one variable
  named `str_endLine` (c-string). For printf, do
  `printf(outFILE, "text%s", str_endLine);`. To add
  line ending to variable do.

```
buffStr[indexSL++] = str_endLine[0];`
if(str_endLine[1] != '\0')
   buffStr[indexSL++] = str_endLine[1];`
```

You can set the line ending with flags

| Flag              | OS         | ending |
|:------------------|------------|--------|
| `-DWINDOWS`       | windows    | `\r\n` |
| `-DLINUX`         | linux      | `\n`   |
| `-DOLD_MAC`       | 90s mac    | `\r`   |
| `-DWIN_BREAK`     | any        | `\r\n` |
| `-DLINUX_BREAK`   | any        | `\n`   |
| `-DOLD_MAC_BREAK` | do not use | `\r`   |

Table: table of unix compiler flags to set the line break
  style. For windows replace `-D` with `/D`. Also never
  use `-DOLD_MAC` or `-DOLD_MAC_BREAK`. These are ancient
  (90s) systems that are no longer used.

**Example**

```
printf("hellow world!%s", str_endLine);
  /*will change line break based on compiler flags*/
```
