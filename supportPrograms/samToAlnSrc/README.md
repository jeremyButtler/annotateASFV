# Use:

Convert sam files to aligned files using the expanded
  cigar format.

The snp output is off when the cigar is not an eqx cigar
  (has 'M' instead of '=' and 'X'). You can switch to
  an eqx cigar with cigToEqxCig.

# License:

Is dual licensed under public domain or MIT. The primary
  license is public domain. However, if public domain
  is not acceptable, will not work, or not wanted, then
  this is under the MIT license.

# Build:

```
make -f mkfile.os # (.linux for bsd/linux; .mac for mac)
make -f mkfile.os install
```

# Use:

You will need to provide the mapped reference file and a
  sam file. Alignments without sequences are ignored. You
  can get the help message with `samToaln -h`.

```
samToAln -ref ref.fasta -sam alignment.sam > out.sam
```

I am using an expanded cigar format becuase I find it
  easier to find the =, X, D, and I symbols to be more
  readable (spot differences faster) then '|' or space
  used by blast and Emboss.

# Example output

I would not recommend using this on reads. It is more for
  the consensus output from tbCon. This is what I had
  on hand.

```
#PG	ID:minimap2	PN:minimap2	VN:2.26-r1175	CL:minimap2 --eqx -a ../freezeTBScripts/TB-NC000962.fa /dev/fd/63
#PG	ID:samToAln	PN:samToAln	VN:2024-07-31	CL:samToAln -ref ../freezeTBFiles/NC000962.fa -sam /home/boom/BC02.sam -out del.sam

##########################################################
# qry: 7b2aa3b6-69a1-443b-a78d-a0d81892abb1
#    aligned: 121 to 341
# ref: NC_000962.3
#    aligned: 3122585 to 3122817
#
# matches: 206
# snps: 7
# insertions: 8
# deletions: 20
# masked: 120
# anonymous: 0
# edit distance: 35
# score: 280
# mean q-score: 20.730206
# median q-score: 20.730206
#
# Eqx legend:
#    = is match
#    X is snp
#    I is insertion
#    D is deletion
#    S is soft mask
##########################################################

ref:   3122585 GTTTCCGTCCCCTCTCGGGGTTTTGGGTCTGACGACCT-CGGA
qry:       121 GTTTCCGTCCCCTC---GGGTTTTGGGTCTGAC-A--TGCGAA
cig:           ==============DDD================D=DD=I==X=

ref:   3122627 CAGCATCTCCCCGGGCGGGCAGCA-----GATATCCCA-TGTT
qry:       158 -------TCC----GCTGTCAGCACATGGGAT-TCCGAGTGTT
cig:           DDDDDDD===DDDD==X=X=====IIIII===D===X=I====

ref:   3122664 TCCGTCCCCTCTCGGGGTTTTGGGTCTGACGACCGACCCGTGG
qry:       189 TCCGTCCCCTCTCGGGGTTTTGGGTCTGACGGTCGACCCGTGG
cig:           ===============================XX==========

ref:   3122707 CCGCCAG-GTTGCCGCCGCCGTTGCTCACCTGGTTTCCGTCCC
qry:       232 CCGCCGGAGTTGCCGCCGCCGTTGCTCACCTGGTTTCCGTCCC
cig:           =====X=I===================================

ref:   3122749 CTCTCGGGGTTTTGGGTCTGACGACCCGGAAGTCAACTAGAGC
qry:       275 CTCTCGGGGTTTTGGGTCTGACGACCCGGAAGTCAACTAGAGC
cig:           ===========================================

ref:   3122792 GGGTGTCGAACGCTGCCCGGTTTCCG
qry:       318 -GGTGTCGAACGCTG-CCGGTTTCCG
cig:           D==============D==========
```
