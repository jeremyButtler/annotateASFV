CC=cc
CFLAGS=\
   -DBYTEMATRIX\
   -DINSDELSNP

python:
	CFLAGS="$(CFLAGS)" CC=$(CC) python3 setup.py install;
	rm -r build;
	rm -r dist;
	rm -r alnSeq.*egg*;

pythonlocal:
	CFLAGS="$(CFLAGS)" CC=$(CC) python3 setup.py build;
	pip3 install --user .;
	rm -r build;
	rm -r alnSeq.*egg*;
	#rm alnSeq.*.so;

debug:
	CFLAGS="$(CFLAGS) -O0" CC=$(CC) python3 setup-debug.py build;
	pip install --user .;
	rm -r build;
	rm -r alnSeq.*egg*;
	egdb -ex "r alnSeqPython.py" -ex "b pyAlnSeq.c:360" --args python3;   
	#gdb -ex "r alnSeqPython.py" -ex "b pyAlnSeq.c:360" --args python3;   

