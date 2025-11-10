/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' seqST SOF: Start Of File
'   - holds gz functions for seqST
'   - the reason this is separate is because my inflate
'     libraries adds around 18kb to all programs, even
'     ones that do not use it
'   o header:
'     - guards and defined variables
'   o fun01: get_gzSeqST
'     - gets a sequence from a fastx/fastx.gz file
'   o license:
'     - licensing for this code (public domain / mit)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|  - guards and defined variables
\-------------------------------------------------------*/

#ifndef GZ_SEQUENCE_STRUCTURE_H
#define GZ_SEQUENCE_STRUCTURE_H

struct seqST;
struct file_inflate;

#define def_noType_gzSeqST 0
#define def_faType_gzSeqST 1
#define def_fqType_gzSeqST 2
#define def_gzType_gzSeqST 4

/*-------------------------------------------------------\
| Fun01: get_gzSeqST
|   - gets a sequence from a fastx/fastx.gz file
| Input:
|   - fileSTPtr:
|     o file_inflate struct to use in uncompression
|   - typeSCPtr:
|     o is set to or has the correct file type
|   - seqSTPtr:
|     o seqST struct pointer to add new sequence to
|   - inFILE:
|     o FILE pointer to file to add (new file)
|     o 0 if already called get_gzSeqST for this file
|     o must be opened in binary format
|       * fopen(<name>, "rb");
|     o you should set inFILE to 0/null after calling,
|       because freeStack_file_inflate 
|       freeHeap_file_inflate will close inFILE for you
| Output:
|   - Modifies:
|     o if have inFILE, sets fileSTPtr be blanked and then
|       have inFILE
|     o if have inFILE, typeSCPtr to have file type 
|       o def_noType_gzSeqST if end of file
|       o def_faType_gzSeqST if fasta file
|       o def_faType_gzSeqST | def_gzType_gzSeqST if
|         fasta.gz file
|       o def_fqType_gzSeqST if fastq file
|       o def_fqType_gzSeqST | def_gzType_gzSeqST if
|         fastq.gz file
|     o fileSTPtr to be on next read/sequence
|     o seqSTPtr to have next read/sequence
|   - Returns:
|     o 0 for no errors
|     o def_EOF_seqST if at end of file
|     o def_memErr_seqST for memory errors
|     o def_badLine_seqST | def_fileErr_seqST for invalid
|       file entries
|     o def_fileErr_seqST for file errors (including gzip
|       checks)
\-------------------------------------------------------*/
signed char
get_gzSeqST(
   struct file_inflate *fileSTPtr, /*has file to extract*/
   signed char *typeSCPtr,         /*gets/has file type*/
   struct seqST *seqSTPtr,         /*gets sequence*/
   void *inFILE
);

#endif

/*=======================================================\
: License:
: 
: This code is under the unlicense (public domain).
:   However, for cases were the public domain is not
:   suitable, such as countries that do not respect the
:   public domain or were working with the public domain
:   is inconvient / not possible, this code is under the
:   MIT license.
: 
: Public domain:
: 
: This is free and unencumbered software released into the
:   public domain.
: 
: Anyone is free to copy, modify, publish, use, compile,
:   sell, or distribute this software, either in source
:   code form or as a compiled binary, for any purpose,
:   commercial or non-commercial, and by any means.
: 
: In jurisdictions that recognize copyright laws, the
:   author or authors of this software dedicate any and
:   all copyright interest in the software to the public
:   domain. We make this dedication for the benefit of the
:   public at large and to the detriment of our heirs and
:   successors. We intend this dedication to be an overt
:   act of relinquishment in perpetuity of all present and
:   future rights to this software under copyright law.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICSLAR PURPOSE AND NONINFRINGEMENT.  IN NO
:   EVENT SHALL THE AUTHORS BE LIABLE FOR ANY CLAIM,
:   DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
:   CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR
:   IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
:   DEALINGS IN THE SOFTWARE.
: 
: For more information, please refer to
:   <https://unlicense.org>
: 
: MIT License:
: 
: Copyright (c) 2024 jeremyButtler
: 
: Permission is hereby granted, free of charge, to any
:   person obtaining a copy of this software and
:   associated documentation files (the "Software"), to
:   deal in the Software without restriction, including
:   without limitation the rights to use, copy, modify,
:   merge, publish, distribute, sublicense, and/or sell
:   copies of the Software, and to permit persons to whom
:   the Software is furnished to do so, subject to the
:   following conditions:
: 
: The above copyright notice and this permission notice
:   shall be included in all copies or substantial
:   portions of the Software.
: 
: THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF
:   ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT
:   LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
:   FOR A PARTICSLAR PURPOSE AND NONINFRINGEMENT. IN NO
:   EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE
:   FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN
:   AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
:   FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE
:   USE OR OTHER DEALINGS IN THE SOFTWARE.
\=======================================================*/
