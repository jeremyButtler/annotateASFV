# Use:

Reverse complements sequences in a fasta file. Not very
  great, but again is here as a individual tool. Many
  repos, such as Emboss have this as well.

# License:

Primary licesenes is the Unlicense. However, if that will
  not work or is not desired for any reason, then this
  defaults to the MIT license.

# Install:

## Unix:

```
if [ ! -d /usr/local/bin ];
then
   sudo mkdir -p /usr/local/bin;
fi

cd ~/Downloads;
git clone https://github.com/jeremyButtler/bioTools;
cd bioTools/revCmpSrc;
make -f mkfile.unix;
suod make -f mkfile.unix install;
```

## Windows:

Not tested, but after downloading bioTools to downloads
  and unzipping. You will need to open a developer
  terminal (it is in visual studio build tools).

```
cd "%HOMEPATH%\Downloads\bioTools\revCmpSrc"
nmake /F mkfile.win
```

Then install revCmp.exe to your target location. Remeber,
  this is a command line program, so you will have to call
  it from command line.

# Using:

Do `revCmp -h` to get the help message.

Do `revCmp -out out.fa file_1.fa` to save the reverse
  complement of file_1.fa to out.fa.

Do `revCmp -out out.fa file_1.fa file_2.fa file_3.fa ...`
  to save the reverse complement of all input fasta files
  to out.fa.

You can reomve `-out out.fa` to output to stdout
  (terminal).
