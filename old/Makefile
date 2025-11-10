all:
	make -C alnSeqSrc && mv alnSeqSrc/alnSeq alnSeq;
	make -C alnSeqSrc CFLAGS="-DDELINSSNP" && mv alnSeqSrc/alnSeq alnSeqMem;
	make -C alnSeqSrc/20230908-alnSeq/ && mv alnSeqSrc/20230908-alnSeq/alnSeq alnSeqScan;
