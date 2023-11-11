all: ssw bio-alignment sequence_align biopython alnSeq

CC=gcc
CXX=g++

ssw:
	git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library || printf "ssw already downloaded\n";
	cd Complete-Striped-Smith-Waterman-Library/src || printf "error happened when downloading ssw\n";
	make CC=$(CC) CXX=$(CXX) -C Complete-Striped-Smith-Waterman-Library/src || gmake CC=$(CC) CXX=$(CXX) -C Complete-Striped-Smith-Waterman-Library/src;

bio-alignment:
	git clone https://github.com/masyagin1998/bio-alignment;
	make CC=$(CC) -C bio-alignment || gmake CC=$(CC) -C bio-alignment;

seq-align:
	git clone --recursive https://github.com/noporpoise/seq-align;
	make CC=$(CC) -C seq-align || gmake CC=$(CC) -C seq-align;

biopython:
	pip install biopython --update;

sequence_align:
	pip install sequence_align;

alnSeq:
	make -C ../ benchmark;
	mv ../alnSeqTwoBit ./;
	mv ../alnSeqByte ./;
	mv ../alnSeqMid ./;
	mv ../alnSeqFast ./;

clean:
	rm -r -f bio-alignment || printf "";
	rm -r -f Complete-Striped-Smith-Waterman-Library || printf "";
	rm -r -f seq-align || printf "";
	rm alnSeqOFast || printf "";
	rm	alnSeqO3 || printf "";
	rm	alnSeqO2 || printf "";
	rm	alnSeqO0 || printf "";
	rm	alnSeqHirschByte || printf "";
	rm	alnSeqHirschTwoBit || printf "";
	rm alnSeqByte || printf "";
	rm alnSeqTwoBit || printf "";
	rm alnSeqFast || printf "";
	rm alnSeqMid || printf "";
	pip remove sequence_align || printf "";

# Requires emboss be installed separately
run:
	bash 20231022-benchmarking/20231022-benchAln.sh;
	Rscript 20231022-benchmarking/20231022-alnSeq-benchmarking.r;
