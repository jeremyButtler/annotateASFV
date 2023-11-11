from alnSeq import alnSeqWater
import sys

#from cbreakpoint import cbreakpoint
#cbreakpoint(breakpoint_id=1)

# Get files and open files
refFileStr = "../analysis/genomes/Small-ref.fasta";
qryFileStr = "../analysis/genomes/Mid-query.fasta";
refFILE = open(refFileStr, "r");
qryFILE = open(qryFileStr, "r");

refHeadStr = "";
qryHeadStr = "";
refSeqStr = "";
qrySeqStr = "";
tmpStr = "";

refHeadStr = refFILE.readline();  # Move past the header
refHeadStr = refHeadStr.replace("\n", "");
refSeqStr = ''.join(refFILE.readlines()).replace("\n", "");
   # readlines reads in all remaining lines in the file
   # ''.join converts the list to a string
      # the '' is needed to fool python into thinking it
      # is a string
   # .replace("\n", "") is to remove new lines

#re.sub("\n", "", refSeqStr); # for regular expressions

qryHeadStr = qryFILE.readline();  # move past the header
qryHeadStr = refHeadStr.replace("\n", "");
qrySeqStr = ''.join(qryFILE.readlines()).replace("\n", "");
   # readlines reads in all remaining lines in the file
   # ''.join converts the list to a string
      # the '' is needed to fool python into thinking it
      # is a string
   # .replace("\n", "") is to remove new lines

alnList = alnSeqWater(gapOpen=-10, fullAln=False, ref=refSeqStr, query = qrySeqStr);
#alnList = alnSeqWater(gapOpen=-10, ref=refSeqStr, query = qrySeqStr);

print(refHeadStr);
print(alnList[0]);

print(qryHeadStr);
print(alnList[1]);

print(alnList[2]);
