# Goal:

describe my expericne and rules for OS portable coding.

# Limitations and my observations

It looks like most languages work well, so long as you
  are on a major OS (most Unix's (Mac, Linux, Bsd, ...)
  and Windows). However, the problems arise when you
  start trying to go to odder OS's (like Plan9). In
  these cases it looks like C89 (C using the 1989
  standard) seems to be present in the high majority
  of cases.

Not a great table, but list of langauges I think are
  on the OS's. 

| language | Windows | Linux | BSD | Mac | Plan 9 |
|:---------|:--------|:------|:----|:----|:-------|
| c89      | cl.exe  | cc    | cc  | cc  | default|
| c99      | no      | cc    | cc  | cc  | no (?) |
| c11      | cl.exe  | cc    | cc  | cc  | no (?) |
| c++      | yes     | yes   | yes | yes | maybe  |
| rust     | yes     | yes   | yes | yes | no (?) |
| go       | yes     | yes   | yes | yes | yes    |
| lisp     | yes     | yes   | yes | yes | sheme9 |
| python 2 | yes     | yes   | yes | yes | yes    |
| python 3 | yes     | yes   | yes | yes | no (?) |
| java     | yes     | yes   | yes | yes | ?      |
| R        | yes     | yes   | yes | yes | no     |
| tcl      | yes     | yes   | yes | yes | yes    |
| tcltk    | yes     | yes   | yes | yes | no     |
| assembler| yes     | yes   | yes | yes | yes    |
| haskell  | yes     | yes   | yes | yes | yes    |

Table: list of languages I think are on variaous OSs
  I know I missed most of them. For C is specifiy if
  you have to install. On many Linux's you will need
  to install cc. No (?) means suspect no, but it may
  be supported. ? means no idea. cl.exe must be
  installed on Windows.

Also, I am assuming the OS can do dynamic memory
  allocaiton and is at least 32 bit (most biology
  programs would struggle on 16 bit). Also, embeded OS's
  without C standard libraries are a differnt game (no
  experince).

Extra features, such as file detection (beyond opening
  user provided file names) are not always the same.
  This also includes threading, SIMD (not always
  supported), GPU acceleration (not always supported),
  and GUIs (never the same).

Again for major OS's you can often really on some of these
  features (openBSD does not support AVX2 and AVX512).
  There are plenty of GUI libraries that are portable
  across multiple OS's. However libraries may vary between
  OS's (not much experience here).

Often the best solution I found in these cases are to
  always have a scalar or command line (for GUIs) fall
  back. However, since characters like carriage return
  (`\r`) can have different functions on different OS's,
  fancy termianl TUI's may not be the most realiable. To
  be honest, the best bet is to really on something that
  does not try to overwrite old lines (like ed).

| OS      | carriage return output                     |
|:--------|:-------------------------------------------|
| Windows | No idea, never done                        |
| Mac     | likley same as Linux                       |
| Old Mac | new line                                   |
| RISC    | no idea (is it even still used?)           |
| Linux   | go to start of line (allows overwrite old) |
| Plan9   | print carriage return character (painfull) |

Table: different printed values for a carraige return.
  Yes plan9 litterally prints a special carriage return
  symbol. Old Mac is pre-Unix days (pre-2000).

| OS      | line break | OS status      |
|:--------|:-----------|:---------------|
| Windows | `\r\n`     | common         |
| Mac     | `\n`       | common         |
| Old Mac | `\r`       | not used       |
| RISC    | `\n\r`     | think not used |
| Linux   | `\n`       | common         |
| Plan9   | `\n`       | hobby OS       |

Table: character(s) for each OS's line break. By hobby
  OS I mean research/hobby, so unlikely to be seen.

To avoid new line and carraige return issues I have
  started to look for all cases. However, for printing,
  Windows will treat `\r\n` as two lines. So, only print
  a `\n`.

```
/*get off line endings*/
if(bufferStr[index] == '\n')
{ /*If: new line line break*/
   ++index;

   if(bufferStr[index] == '\r')
      ++index;
} /*If: new line line break*/

else if(bufferStr[index] == '\r')
} /*Else If: carriage return line break*/
   ++index;

   if(bufferStr[index] == '\n')
      ++index;
} /*Else If: carriage return line break*/
```

# Build systems

## CMake

Sounds like this is what cmake was designed for. Maybe
  the best bet for Unix (Linux/Mac/Bsd) and Windows. No
  idea if can target plan9.

## Make

I use a script to build my make files. That way when I
  need to upadate make files after adding dependencies,
  I can just update the script and then rebuild
  everything.

### Make on Unix

Is not realiable and even on Unix's vary a lot. For
  example Mac and Linux use gnu-make, while BSDs use
  bsd-make. On both OS's the native make is called
  make. So, when building a portable make file you
  must limit yourself to orignall make. The major
  parts are no if statments and no subshell
  commands (`VARIABLE=$(shell command)`).

Of course BSDs can get gnu-make (`gmake`), but gmake
  then reverts back to bsd-make when you call make within 
  make `make -F /path/to/mkfile -C /path/to/build`. The
  best solution I have found for this is to either write a
  make file that works with both systems or to include a
  make variable `MAKE=make` and then call make by the
  variable `$(MAKE) MAKE=$(MAKE) -F /path/to/mkfile -C /path/to/build`.

```
LD=cc
CC=cc
CFLAGS=-DNONE # so user can modify
coreCFLAGS= -O3 -c -std=c89 <you input ...>
   # do no want user modifying
NAME=program
PREFIX=/usr/local/bin
O=o.unix

objFiles= \
   file1.$O \
   file2.$O \
   mainFile.$O;

program: mainFile.$O
	$(LD) -o $(NAME) $(objFiles)

mainFile.$O: mainFile.c file1.$O file2.$O
	$(CC) -o mainFile.$O $(CFLAGS) $(coreCFLAGS) mainFile.c

file1.$O: file1.c file1.h
	$(CC) -o file1.$O $(CFLAGS) $(coreCFLAGS) file1.c

file2.$O: file2.c file2.h
	$(CC) -o file2.$O $(CFLAGS) $(coreCFLAGS) file2.c

install
	if [ ! -d "$(PREFIX)" ];
   then
      mkdir "$(PREFIX)";
	fi # deal with mack not having /usr/local/bin
   mv "$(NAME)" "$(PREFIX)";
   chmod a+x "$(PREFIX)/$(NAME)";
```

## Make on Windows

For Windows, it gets even worse, because windows `nmake`
  has its own syntax for if statements and the complier
  c flags are different. I typically will create a
  separate make file for windows. Also, for really easiy
  installs you should make a `.bat` script to do the
  install or figure out how to make an `.msi. file.

Also Windows `nmake` does not have `-C`, so you need to
  use `cd /path/to/build` to move to your build directory
  and then `nmake /F /path/to/mkfile`.

```
LD=link.exe
CC=cl.exe
NAME=program.exe
CFLAGS=/DNONE # so user can modify
coreCFLAGS= /c /O2 /Ot /Za /Tc
   # do no want user modifying
O=o.win

objFiles= \
   file1.$O \
   file2.$O \
   mainFile.$O;

program: mainFile.$O
	$(LD) -o $(NAME) $(objFiles)

mainFile.$O: mainFile.c file1.$O file2.$O
	$(CC) -o mainFile.$O $(CFLAGS) $(coreCFLAGS) mainFile.c

file1.$O: file1.c file1.h
	$(CC) -o file1.$O $(CFLAGS) $(coreCFLAGS) file1.c

file2.$O: file2.c file2.h
	$(CC) -o file2.$O $(CFLAGS) $(coreCFLAGS) file2.c
```

Another problem is that windows does not have the
  compliers storted in the PATH. So, you will need
  to either have the user fire a developer terminal or,
  have a batch (bat) script activate it for you. Here is
  an example of a bat script I used to acitivate the
  developer enviroment for the visual stuido 2022
  build tools.

```
@echo off

:: pushes directory move to stack so can undo later
pushd %~dp0

set guiDir="..\programs\guiFreezeTBSrc"

:: check if in a development shell (have cl.exe) or
:: try to actiavate the development shell settings

:: from stack overflow; find exectuable in path
::   if cl.exe is in path, then in a developer enverioment

for %%X in (cl.exe) do (set FOUND=%%~$PATH:X)

:: this is the if statment to acivate the developer
::   enviroment. You can add in extra if's to add extra
::   paths (maybe?)
if defined FOUND (
   echo "in developer envrioment"
) else (
   if exist "%programFiles(x86)%\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat" (
      CALL "%programFiles(x86)%\Microsoft Visual Studio\2022\BuildTools\VC\Auxiliary\Build\vcvars64.bat"
) else (goto :error))

:: the CALL command has called bat script operate in the
:: current bat scripts terminal, otherwise the other
:: script fires a separate terminal

:: compile code
nmake /F mkfile.win
nmake /F mkfile.win clean

:: see if can install in program files (this file does
:: not exit)
copy /y nul "%programFiles%\2025-02-11--test-file.txt"

:: install program in app data (global) or program
if exist "%programFiles%\2025-02-11--test-file.txt" (
   :: remove test file I did to check the path
   del "%programFiles%\2025-02-11--test-file.txt"

   move program.exe "%programFiles%" /E /I
) else (
   move program.exe "%localAppData%" /E /I
)
```

### Make on Plan9

Plan9's `mk`, like everything else in Plan9 is odd. The
  compiler and linker are auto-detected by the system.
  Also, variables do not use `()`. Otherwise, It is very
  similar to make otherwise. You could likley use `ape`
  (Plan9's old c89 compiler) instead of the native
  compilers.

```
</$objtype/mkfile
   # defines as make file type and gets varaibles

CFLAGS=-DNONE # so user can modify
coreCFLAGS=-DPLAN9 # do no want user modifying
NAME=program
PREFIX=$home/bin # or /bin (no idea?)
O=o.plan9

objFiles= \
   file1.$O \
   file2.$O \
   mainFile.$O;

program: mainFile.$O
	$LD -o $NAME $objFiles

mainFile.$O: mainFile.c file1.$O file2.$O
	$CC -o mainFile.$O $CFLAGS $coreCFLAGS mainFile.c

file1.$O: file1.c file1.h
	$CC -o file1.$O $CFLAGS $coreCFLAGS file1.c

file2.$O: file2.c file2.h
	$CC -o file2.$O $CFLAGS $coreCFLAGS file2.c

install:
   mv $NAME $PREFIX;
   chmod a+x $PREFIX/$NAME; # no idea if needed
```

## my own make

Hopping to get up one of these days. Likley cmake will
  be better.

# C89 differences

To get C protable with the native compilers on Plan9 I
  have to use a slightly stricter standarded then just
  C89. However, I should make some notes about the c89
  standerd. Mainly becuase with `gcc` it is easy to use
  gnu extensions without realizing it.

- Comments in c89 are only block comments (`/**/`), so
  `\\` is not allowed

- All variables must be declared at the start of the
  function. 
  - this includes variables in loops,
    so `for(int x = 0; x < 10; ++x){...}` is not allowed

- Macros must be on a single line
  - gcc will not enfource this with `-std=c89` only
  - Windows also requires this

Do not do:

```
#define brillentMacro(x) \
( \
   (x) + 10 \
)

```

Instead do:

```
#define brillentMacro(x) ((x) + 10)
```

Or if only use a few times, directly inline the macro.

- Arrays on stack must be defined by constants
  - so `int x = 10; char string[x];` is not allowed
  - instead `int x = 10; char string[10];`
  - instead `#define x_fun10 10 char string[x_fun10];`

- `#define variable` are effectivly global scope, so I
  can not use two defines in two separate functions with
  the same name (plan9 enfources this)

Do not do this:

```
fun01(){
   #define y 10
   return y;
}

fun02(){
   #define y 10
   return y;
}
```

Instead do this:

```
fun01(){
   #define y_fun01 10
   return y_fun01;
}

fun02(){
   #define y_fun02 10
   return y_fun02;
}
```

- You can only use the c89 libraries. The ones I have used
  are `<stdlib.h>`, `<sdtio.h>`, `<time.h>`, `<string.h>`.
  There are a few more, so check the c-89 draft standard.
  - You can include libraries from others, but make sure
    that they are at least as protable as your code aims
    to be.
  - Avoid non-standard system libaries, or have
    fallbacks/alternatives for your target OS's

There is probably more, but this is what stands out to me.

# OS quirks

## Mac:

Most stuff works well between linux and Mac.

## Windows:

The linux `popen` library is under a different name. As a
  rule do not use `popen`. Do not confuse `popen` with
  `fopen`. `popen` is to open a new pipe (not stdout),
  `fopen` is for files.

The only bug I know on windows is never use `fgets` with
  `fseek` or `ftell`.

Also, in Linux fread functions like a binrary reader,
   however, in Windows fread does not. Always use
   `fopen(<file>, "rb");` when you want to read a file
   byte by byte (likely keeps `\r\n` line endings) instead
   of getting c-strings.

DO NOT DO:

```
temporaryPointer = fgets(buffer, 100, inFILE);
position = ftell(inFILE);
```

```
temporaryPointer = fgets(buffer, 100, inFILE);
position = fseek(inFILE, 25, SEEK_SET);
```

Instead (might have sizeof position wrong):

```
readIn = fread(buffer, sizeof(char) , 100, inFILE);
position = ftell(inFILE);
```

```
readIn = fread(buffer, sizeof(char), 100, inFILE);
position = fseek(inFILE, 25, SEEK_SET);
```

Other then bugs, you should be aware that windows
  uses `\r\n` for line endings. So, when reading data
  look for both `\r` and `\n`. Also account for the
  possibility of them bing together.

## Plan9

Plan9 (with native compilers [not ape]) makes the game a
  harder, because it forces a coding style on you.
  However, most of it is just modifying your coding style
  a bit.

No idea if even has a `popen` equivlent. Do not confuse
  `popen` with `fopen`. `popen` is to open a new pipe (not
  stdout), `fopen` is for files.

### Plan9 library includes (only .c files)

The first and hardest is to not use `#include` in
  your `.h` files (`.c` files only). You can get
  around structs by forward declaring them in your
  `.h` files (`struct structName`). For stdio.h's FILE
  structure, you will have to pass `void *inFILE` to your
  functions. Then use `(FILE *) inFILE)` when using the
  file.

For foward declarations; all structures must be declarated
  as pointers, since the C compiler does not know the
  structure size until runtime.

example of forward declaration in samEntry.h file (setting up
  struture):

```
typedef struct samEntry{
   signed char *sequence;
   signed int sequenceLength;
   signed int sequenceSize;
}samEntry;

/*initialization functions go here*/
```

example of forward declaration in bioTools.h file:

```
struct samEntry;

typedef struct bioTools{
   <variables>

   struct *samEntry *samPointer;
     /*do not do: struct samEntry samStack;*/
}bioTools;

void
printSam_bioTools(
   struct bioTools,
   void *samFILE
);

signed char
addSequenceToSamEntry_bioTools(
   signed char *sequence,
   struct *samEntry
);
```

example of using foward delare in bioTools.c file.

```
#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
   #include <string.h>
#endif

#include "bioTools.h"
#include <stdio.h>
#include "samEntry.h"

void
printSam_bioTools(
   struct bioTools,
   void *samFILE
){
   fprintf(
      (FILE *) samFILE,
      "%s\n",
      bioTools->samPointer->sequence
   )
}

signed char
addSequenceToSamEntry_bioTools(
   signed char *sequence,
   struct *samEntry
){
   signed int length = strlen(sequence);

   if(length >= samEntry->sequenceSize)
   {
      free(samEntry->sequence);
      samEntry->sequence = 0;
      samEntry->sequence =
          malloc((length + 1) * sizeof(signed char))

      if(! samEntry->sequence)
         return 1; /*memory error*/

      samEntry->sequenceSize = length;
   }

   strcpy(
      samEntry->sequence,
      sequence
   );

   samEntry->sequenceLength = length;
}
```

### Plan9 Pragmas, limited

For pragmas, you can only use `#ifdef`, `#ifndef`,
  `#else`, and `#endif`. This means that you will be doing
  a lot of nesting. You can sometimes reduce it down if
  you konw only one case will be fired.

nested cases:

```
#ifdef CASE_ONE
   doe something for case one
#else
   #ifdef CASE_TWO
         doe something for case two
      #else
         doe something for case three
      #endif
   #endif
#endif

```

If I know only one case will ever be provided, I can
  avoid nesting.

```
#ifdef CASE_ONE
   doe something for case one
#endif

#ifdef CASE_TWO
   doe something for case two
#endif

#ifdef CASE_THREE
   doe something for case three
#endif
```

### Plan9 library names

Plan9 has also decided to merge most of the standard
  libraries into libc and also uses a u.h library. This is
  a bit of pain, because it means I have to have a special
  header for it. The only expection I know of is stdio.h,
  which must come after libc (which must come after u.h).

So as a rule, if it a non-stdio.h standard library, then
  assume it is in libc.h.

```
#ifdef PLAN9
   #include <u.h> /*always comes beforre libc*/
   #include <libc.h> /*always comes beforre libc*/
#else
   #include <stdlib.h>
   #include <string.h>
   #include <time.h>
#endif

#include <stdio.h>
```

To be honest my main problem with this is that it is not
  standard, otherwise it would be great.

## Plan9 long is long long

In plan9 they authors decided to keep longs as 32 bits, so
  inorder to get 64 bit longs you will need to do
  `long long`. My trick here is to declare a plan9 long
  variable and use that when I want the largest value
  possible on the system.

```
#ifdef PLAN9_64
   typedef unsinged long long ul_filename;
   typedef singed long long sl_filename;
#endif
```

I use this trick rarely, but you can see it in ulCp and a
  few other of my libraries. Just keep in mind you can not
  use these in the .h file, unless it is the .h file that
  defines these. So, each file gets its own plan9 long
  when I needed them.

# Plan9 return not exit

For Plan9 the exit function returns `void *`. So, using
  `exit(1)` will mess you up. A solution is instead to
  use `return(1)` and `return(0)`.

```
int
main(
   int numArgs,  /*always is int to avoid warings*/
   char *args[]  /*always is char to avoid warings*/
){
   /*add code*/

   if(error)
      return error;
   else
      return 0;
}
```

## CPU quirks

Here I am very weak. I only have experince on little endin
  systems.

ARM CPUs treat characters as unsigned, while x86 CPUs
  treat characters as signed. So, for characters (at
  least) always declare your sign (
  ex: `signed char *buffer` or `unsigned char *buffer`).

### Some basic tricks to get values

The best reasource for bit hacking is the stanford bit
   hacks guide
   [https://graphics.stanford.edu/~seander/bithacks.html](
    https://graphics.stanford.edu/~seander/bithacks.html).

You will not have access to the stdint library and mabye
  even the limits libraries. However, there are some
  tricks you can do (stick this in your plan9 64bit .h
  file).
 
- maximum value of unsigned datatype
  - ex: `#define utypeMax ( (usigned type) -1 )`.
- maximum value of signed datatype
  - ex: `#define stypeMax ( (signed type) (((usigned type) -1) >> 1) )`.
- number bytes in type
  - ex: `#define typeBytes (sizeof(type))`
- number bits in unsigned type
  - ex: `#define utypeBits (sizeof(type) << 3)`
- number bits in signed type
  - ex: `#define stypeBits ( (sizeof(type) << 3)`>> 1)
- truncate pattern to type to target bit size
  - ex: `#define pattern (unsigned long) 0xf1f1f1f1f1f1f1f1`
  - at this point pattern is reduce to the size of an
    unsigned long on the system
- negative bit
  - ex: `#define typeNegBit ((unsigned type) 1 << (utypeBits - 1))`
  - at this point pattern is reduce to the size of an
    unsigned long on the system

You can then use these pre-defied variables in loops to
  have a variable opertions. Since the loops have
  hardcoded lengths, -O3 will remove. Also, you can
  often use the pre-defined variable in your bit
  operations.

```
for(
   bitUL = 0;
   bitUL < (sizeof(unsigned long) << 3);
   ++bitUL
){
   /*perorm operation on bit*/
}
```

### SIMD, is missing

You will need a fallback for SIMD, since not all OS's
  support it. However, in some limited cases you can use
  longs to copy strings, for addition/subtraction, and
  maybe some basic bitwise (non-shift) operations. Keep in
  mind, for additon/subtraction, an overflow will mess you
  up. So, make sure there is no chance of overflow.

Also, this trick only works on CPUs/OSs that can handle
  odd byte alignments for ints/longs. So, you should
  always have a byte fallback for those cases (I think
  it is RISC OS?, I know modern ARM CPUs are not an
  issue). I use the `-DNOUL` flag to enable byte copying
  mode.

See ulCp for long string coping.

You will want to be a bit carefull here, since this can
  do intresting this with big ending systems. So, as a
  rule, only shift bits to check values (ex: is null).

# Good habits

Some habits I have learned. These may be good, or may be
  be bad. Some of these most people do, but I had to learn
  these on my own, so I mention them here.

- Always include `#ifndef FILE_H`, and `#define FILE_H` at
  the top of your .h file and `#endif` at the end
  - I make this note because when coding for Plan9, you
    do not need to do this, however, you never know what
    others will do, so keep it safe
- My experience with global variables in .h files is
  that they generate warnings when not used, so what
  I have learned is to:
  - are only include one `static` vairable per .h file
  - the .h file is dedicated to its one static variable or
    array
    - put functions in a separate .c file
  - global strings can not be more then 512 characters
    long
- When declaring variables, always use signed or unsigned
  - for characters this is a must, since ARM CPUs treat
    characters as unsigned, while x86 CPUs treat
    characters as signed.
- Try to reduce the amount of memory allocated in
  functions
  - each function (except free functions) only frees
    the memory it has alloacted, so that only main
    freese everything
  - when memory is reused, pass buffers in and resize
    buffers to avoid many malloc calls
- Make sure all function allocated varaibles are cleaned
  up or returned at the end of a function
  - This reduces memory leaks and often is little in cost
    compared to the function exection time
  - Also allows making errors into
    gotos `if(error) goto err;`, which then allows me to
    not worry about what I have to free when I have
    errors.
- after debugging (or during), always test code with
  valgrind (sometimes I forget)
- always have at least two make files, one for debugging
  and one for the user
  - I always use `-O0`, `-g`, and `-Werror` in the
    dubugging make file.
  - I never use debugging symbols (`-g`/`-ggdb`) in the
    users make file. The debugging symbols can make a
    400kb program into a 900kb program.
  - most users know this, but I always use `-std=c89`,
    `-Wall`, `-Wextra`, and `-Wpendantic` in both my
    debugging and non-debugging make files
  - treat warnings seriously, unless there is no way
    around or you are on Windows or Plan9, were warnings
    are either to sensitive (Plan9 with extra warnings) or
    about the standard libraries (Windows with extra
    warnings). Tread both Windows and plan9 on a case by
    case bases.
  - ignore openBDS's warning about `strcpy`. Yes it is
    important, but there is no good replacement, including
    openBSD's suggestion on all OSs.
    - my solution was to make my own function that openBSD
      has never seen. Just as problematic, but the
      compiler has never seen it.
      - What else do you expect people to do when you
        complain about the standar libraries?
  - this is a habit I need to learn, remove unused code
    from your .c files or place then in an old .c file.
    they just increase the compile size of your program.
    - the only exception I can think of is if you are
      coding a library of usefull functions.

Here is an example of my clean up.

```
/*not a very good example, because it is to simple*/

getNumbers(x){
   signed int error = 0;
   signed int buffer = malloc(x * sizeof(signed int));

   if(! buffer)
      goto memoryError;

  if(x < 0)
     goto toLowError;

  while(x > 0)
     buffer[x] = x--;

  free(buffer);
  buffer = 0;
  error = 0;
  goto cleanUp;

  memoryError:;
     error = -1;
     goto cleanUp;

  toLowError:;
     error = -2;
     goto cleanUp;

  cleanUp:;
     if(buffer)
        free(buffer);
     buffer = 0;

     return error
}
```

## Weird things I have noticed

One trick I have picked up is that .o and .a files will be
  included in static compile steps, even when you are not
  using the functions that call them. For example:

Say we have library 1:

```
#include "someLib.h"

add(int x, int y){ return x + y;};
use_somLib(int x, int y){ return somLib(x + y);};
```

Even if I never use `use_someLib` in my final program, the
  `someLib` .o files will be included. At one or two
  points I have broken libraries up at logical points to
  prevent this bloat. The two I have done are dealing
  with `.gz` files in gzSeqST instead of seqST and
  splitting `samEntry` into `samEntry` and `samRef`.

I guess what I have learned is if there is a logical
  splits and the extra libraries are rarely used in other
  .c files, then split them into smaller units.

I would avoid using `fileFun`, but that is to critical for
  what I do. Otherwise I could reduce programs the used
  it (most, if not all) by 4kb.

# Summary

OS coding is a pain, mainy because you loose access to
  most libraries. Also, it takes a bit to get used to the
  coding style.

There is a lot I am missing, for example you can
  multithread on Plan9, however, it is differnt and I
  have never dug into it.

GPUs are major OS's only, GUIs are major OS's only,
  unless you willing to build a GUI for additional
  target OSs.

Always make sure you test your code on multipe OS's,
  particullary the differnt ones (plan9/windows/unix).
