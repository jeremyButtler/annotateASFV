# Goal:

Describe how to use the fileFun library to read line
  by line across different OS's.

Files: fileFun.h and fileFun.c

- Dependencies:
  - ulCp
  - endline

# fileFun

## fileFun overview

Non-gz functions for file manipulation.

- [lineCnt\_fileFun](#lineCnt_fileFun)
  - count number of lines in a file
- [getLine\_fileFun](#getLine_fileFun)
  - get a line or maximum number of bytes from a file
- [getFullLine\_fileFun](#getFullLine_fileFun)
  - get a full length line from a file

| Flag           | Use                                   |
|:---------------|:--------------------------------------|
| `-DSLOW_BREAK` | slower, reliable line break detection |

Table: configuartion options. For windows replace `-D`
  with `/D`.

## lineCnt\_fileFun

lineCnt\_fileFun will seek the start of the file after it
  is finished. It will also modify the signed long to hold
  the length of the longest line in the file. This length
  does not include line breaks.

| Input            | Use                                 |
|:-----------------|:------------------------------------|
| `FILE *`         | file to find the number of lines in |
| `signed long *`  | max line length; minus line breaks  |

Table: Input for lineCnt\_fileFun

| Output               |
|:---------------------|
| number lines in file |

Table: output for lineCnt\_fileFun

**Example**

```
FILE *inFILE = fopen("file.txt", "r");
signed long *maxWidthSL = 0;
signed long fileLinesSL =
   lineCnt_fileFun(inFILE, &maxWidthSL);

/*maxWidthSL gets the width of the longest line,
`  fileLinesSL gets the number of lines in the file
`  sets inFILE to be at start
*/
```

## getLine\_fileFun

getLine\_fileFun will return the number of bytes in the
  read in buffer, but it will also modify a signed long
  to have the number of bytes read from the file. These
  values can differ, since getLine\_fileFun will replace
  line breaks (`\n`, `\r\n`, `\r`, and `\n\r`) with the
  OS line breaks it was compiled for. So, the return
  value includes the modified line break ending, while
  the modified signed long has the actual bytes read from
  the file.

| Input            | Use                             |
|:-----------------|:--------------------------------|
| `FILE *`         | file to get next line in        |
| `signed char *`  | buffer to hold output           |
| `signed long`    | maximum number of bytes to grab |
| `signed long *`  | number bytes read from file     |

Table: Input for getLine\_fileFun

| Output | Use                                           |
|:-------|:----------------------------------------------|
| > 0    | number bytes in buffer (can differ from file) |
| 0      | End Of File (EOF)                             |

Table: Output for getLine\_fileFun

**Example**

```
signed char bufferStr[512];
signed long realLengthSL = 0;
signed long numBytesSL = 0;
FILE *inFILE = fopen("file.txt", "r");

numBytesSL =
   getLine_fileFun(inFILE, bufferStr, 521, &realLengthSL);
   /*real length gets the number of bytes fgets read in
   `   numBytesSL gets the number of bytes in bufferStr
   `    - this count excludes the `\r` in `\r\n`
```


## getFullLine\_fileFun

getFullLine\_fileFun will use expand a buffer until a full
  line has been read in. It uses getLine\_fileFun to read
  in line by line.

| Input            | Use                                |
|:-----------------|:-----------------------------------|
| `FILE *`         | file to get next line in           |
| `signed char **` | buffer to hold output (is resized) |
| `signed long *`  | current size of buffer             |
| `signed long *`  | number bytes read from file        |

Table: Input for getFullLine\_fileFun

| Output | Use                                           |
|:-------|:----------------------------------------------|
| > 0    | number bytes in buffer (can differ from file) |
| 0      | End Of File (EOF)                             |
| < 0    | memory error when resizing the buffer         |

Table: Output for getFullLine\_fileFun

**Example**

```
signed char *bufferStr = 0;
signed long bufferSizeSL = 0;

signed long realLengthSL = 0;
signed long numBytesSL = 0;

FILE *inFILE = fopen("file.txt", "r");

numBytesSL =
   getLine_fileFun(
      inFILE,
      &bufferStr,    /*auto allocated*/
      &bufferSizeSL, /*auto filled*/
      &realLengthSL, /*bytes read in; includes '\r'*/
      0              /*start at character 0*/
   );
   /*numBytesSL gets number bytes read, but excludes the
   `  '\r' in "\r\n"
   */

if(numBytesSL)
   /*deal with memory errors*/
```
