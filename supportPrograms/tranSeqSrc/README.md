# Use:

Translate nucleotide sequences into amnio acid sequences.

# Make:

Mac or linux

```
git clone https://github.com/bioTools
cd bioTools/tranSeqSrc
make -f mkfile.unix
sudo make -f mkfile.unix install
```

Windows

```
cd "%HOMEPATH%\bioTools\tranSeqSrc"
nmake /F mkfile.win
```

Copy tranSeq.exe to your install location.

# Call:

You can get the help message with `tranSeq -h`.

You can translate a sequence
   with `tranSeq -fa sequence.fasta > translate.fasta`.
