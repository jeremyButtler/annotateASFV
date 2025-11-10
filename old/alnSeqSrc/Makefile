PREFIX=/usr/local/bin

CC=cc

# These are flags I do not want the user to over write,
# acidentally.
COREFLAGS=\
  -Wall\
  --std=c99\
  -static\
  -O3

# These are here for the user to overwrite
#CFLAGS=-DBLANK
CFLAGS=\
   -DBLANK\
   -DBYTEMATRIX\
   -DINSDELSNP

DEBUGFLAGS=\
   -DNOSEQCNVT\
   -DBYTEMATRIX
   #-DNOGAPOPEN
# Current options: -DNOSEQCNVT, -DHIRSCHTWOBIT

SOURCE=\
  cStrToNumberFun.c \
  seqStruct.c\
  scoresST.c\
  alnSetStruct.c\
  generalAlnFun.c\
  alnMatrixStruct.c\
  alnStruct.c\
  waterman.c\
  needleman.c\
  hirschberg.c\
  memWater.c\
  alnSeq.c
 
# Build findCoInfct
all:
	$(CC) $(COREFLAGS) $(CFLAGS) $(SOURCE) -o alnSeq

python:
	CC=$(CC) make -C pythonPkg/ python;
pythonlocal:
	CC=$(CC) make -C pythonPkg/ pythonlocal;

debug:
	$(CC) -Wall  -static --std=c99 -O0 -ggdb $(DEBUGFLAGS) $(SOURCE) -o alnSeqDebug
	# Used to use -g, but -ggdb provides more info for gdb
	bash debug.sh
	# edit debugCMDs.txt to change the gdb commands

egcc:
	#egcc $(COREFLAGS) $(CFLAGS) $(SOURCE) -o alnSeq
	egcc $(COREFLAGS) $(CFLAGS) $(SOURCE) -o alnSeq
gcc:
	gcc  $(COREFLAGS) $(CFLAGS) $(SOURCE) -o alnSeq
cc:
	cc   $(COREFLAGS) $(CFLAGS) $(SOURCE) -o alnSeq

# These settings are here for quick compiling
fast:
	$(CC) $(COREFLAGS) -DNOGAPOPEN -DBYTEMATRIX -DINSSNPDEL $(CFLAGS) $(SOURCE) -o alnSeqFast

benchmark:
	$(CC) $(COREFLAGS) -DHIRSCHTWOBIT $(SOURCE) -o alnSeqTwoBit
	$(CC) $(COREFLAGS) -DBYTEMATRIX $(SOURCE) -o alnSeqByte
	$(CC) $(COREFLAGS) -DBYTEMATRIX -DINSSNPDEL $(SOURCE) -o alnSeqMid
	$(CC) $(COREFLAGS) -DBYTEMATRIX -DNOGAPOPEN -DINSSNPDEL $(SOURCE) -o alnSeqFast
	CC="$(CC)" make -C pythonPkg/ pythonlocal;

clean:
	rm alnSeqDebug || printf ""; # Only thing to clean up
	rm alnSeqalnSeqHirschByte || printf "";
	rm alnSeqHirschTwoBit || printf "";
	rm alnSeqOFast || printf "";
	rm alnSeqO3 || printf "";
	rm alnSeqO2 || printf "";
	rm alnSeqO0 || printf "";
	rm alnSeq*.core || printf "";
	rm alnSeqByte || printf "";
	rm alnSeqTwoBit || printf "";
	rm alnSeqFast || printf "";
	rm alnSeqMid || printf "";
	# || printf ""; is so it does not error out

install:
	mv alnSeq $(PREFIX) || mv alnSeqFast $(PREFIX) || mv alnSeqMid $(PREFIX) || printf "Unable to install alnSeq at %s\n Change this with make PREFIX=/path/to/install install\n" $(PREFIX) && exit;
	chmod a+x $(PREFIX)/alnSeq*;
