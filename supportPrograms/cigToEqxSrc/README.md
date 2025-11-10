# Use:

Converts normal sam file cigars ("M, I, D, S, H") into eqx
  sam file cigars ("=, X, I, D, S, H").

Really is here for tbCon or in case you forgot to have
  minimap2 output eqx cigars (`--eqx`).

# Install:

Build cigToEqx using `make -f` and the make file (mkfile)
  for your os. The exentstion on the mkfile shows the OS
  it is set up for (mkfile.os).

Make file names:
  - mkfile.unix: any unix (Mac/Linux/BSD) OS
  - mkfile.static: static build (Linux/BSD)
  - mkfile.win: windows (no install command)

```
make -f mkfile.unix
sudo make install
```

# How to use:

To use cigToEqx you will need to provide a sam file and
  the reference that was used to make the sam file. The
  file should only have one reference. You can print the
  help message with `cigToEqx -h`.

```
cigToEqx -sam file.sam -ref ref.fasta > out.sam
```

# Updates

- 2025-01-07:
  - added in multi-ref support
