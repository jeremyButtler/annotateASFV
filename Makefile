all:
	if [ ! -d "$$(pwd)/bin" ]; then mkdir "$$(pwd)/bin";fi;
	make -f mkfile.unix -C supportPrograms;
	make -f mkfile.unix -C supportPrograms \
		PREFIX="$$(pwd)/bin" install;
	cp scripts/annotateASFV.sh bin;
	cp scripts/mapGene.sh bin;
	cp scripts/aaToFeature.awk bin;
	cp scripts/extractSeqByScore.awk bin;
	cp scripts/rmDupFeatures.awk bin;
clean:
	make -f mkfile.unix -C supportPrograms clean;
cleanAll:
	if [ -d "$$(pwd)/bin" ]; then rm -r "$$(pwd)/bin"; fi;
	make -f mkfile.unix -C supportPrograms clean;
