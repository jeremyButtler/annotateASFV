# A good chunk of this is from
# https://github.com/pypa/setuptools/blob/4b64119ab29b1418dc4975ea4212ab38f112c2ab/setuptools/_distutils/ccompiler.py#L964
# Linux:
#  CC=cc python3 setup.py build
#  pip install -e .

from setuptools import setup
from setuptools import Extension

descriptionStr="\
   Neeldeman Wunsch, Waterman Smith, and Hirschberg\
   pairwise alingers written in C for python. See the\
   docstrings for more details.\
";

srcFilesStr = [
      "../alnSeqDefaults.c",
      "../cStrToNumberFun.c",
      "../seqStruct.c", 
      "../scoresST.c",
      "../alnSetStruct.c",
      "../generalAlnFun.c",
      "../alnMatrixStruct.c",
      "../alnStruct.c",
      "../hirschberg.c",
      "../needleman.c",
      "../waterman.c",
      "pyAlnSeq.c"
]; # Source files to compile library
      #" ../memWater.c",

compileFlags = [
   "-O3",
   "-DBYTEMATRIX",
   "-DINSDELSNP",
   "-static",
   "--std=c99"
]
setup(
   name = "alnSeq",
   version = "20231029",
   description = descriptionStr,
   author = "https://github.com/jeremyButtler/alnSeq",
   ext_modules=[
      Extension(
         "alnSeq",
         srcFilesStr,
         extra_compile_args = compileFlags
   )]
);
