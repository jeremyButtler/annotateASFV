import sys # for command line arguments
from sequence_align.pairwise import hirschberg, needleman_wunsch
# https://github.com/kensho-technologies/sequence_align
# sequence_aligne was installed Oct 22, 2023

# Get files and open files
refFileStr = sys.argv[1];
qryFileStr = sys.argv[2];
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

alnRefStr, alnQryStr = needleman_wunsch(
      refSeqStr,
      qrySeqStr,
      match_score=2.0,
      mismatch_score=-1.0,
      indel_score=-2.0,
      gap="-"
);

print(refHeadStr);
print(''.join(alnRefStr));

print(qryHeadStr);
print(''.join(alnQryStr));
