# Use:

Does Waterman alignment of two sequences and returns the
  coordinates and alignment scores of the result. This
  uses less memory than a Waterman, but does not return
  the aligned sequence.

# License:

Dual licensed under the Unlicense or MIT. Pick the one
  that works for you.

# Install:

## Unix (Linux/Mac/BSD):

Only tested on Linux. However, nothing in code to prevent
  portability to Mac.

```
cd ~/Downloads;
git clone https://github.com/jeremybuttler/bioTools;
cd bioTools/memwaterSrc;

if [ ! -d "/usr/local/bin" ];
then
   sudo mkdir -p "/usr/local/bin";
   # this is sometimes missing on MAC
fi

make -f mkfile.unix;
sudo make -f mkfile.unix install;
```

## Windows:

In theory this should work, but is a pain. You will need
  the microsoft C++ build essentials (included with
  visual studio). Then open a developer terminal.

Clone/Download the bioTools repository, unzip and save the
  unzipped files Downloads. Make sure it is named
  bioTools and that you do not have something like
  `bioTools\bioTools`. Then open a developer terminal.

```
cd "%HOMEPATH%\Downloads\bioTools\memwaterSrc"
nmake /F mkfile.win
```

install memwater.exe were you want. Just remember you can
  only run memwater from the terminal.

# Using:

To get the help message do `memwater -h`.

To map two fasta files
  do `memwater -ref ref.fasta -qry query.fasta`.

To map two fastq files
  do `memwater -ref-fq ref.fastq -qry-fq query.fastq`.
