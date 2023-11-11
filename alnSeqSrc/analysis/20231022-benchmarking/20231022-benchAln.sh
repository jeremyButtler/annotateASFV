#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

scriptDirStr="$(dirname "$0")";
smallTestQueryStr="$scriptDirStr/../genomes/Small-query.fasta";
smallTestRefStr="$scriptDirStr/../genomes/Small-ref.fasta";

midTestQueryStr="$scriptDirStr/../genomes/Mid-query.fasta";
midTestRefStr="$scriptDirStr/../genomes/Mid-reference.fasta";

largeTestQueryStr="$scriptDirStr/../genomes/Large-query.fasta";
largeTestRefStr="$scriptDirStr/../genomes/Large-reference.fasta";

ramKillerTestQueryStr="$scriptDirStr/../genomes/RamKiller-query.fasta";
ramKillerTestRefStr="$scriptDirStr/../genomes/RamKiller-ref.fasta";

numRepI=10;
statsFileStr="$scriptDirStr/20231022-alnSeq-bench-stats.tsv"
iCnt=0;

# seq-align
# git clone https://github.com/noporpoise/seq-align
# make

# bio-alignment (has hirschberg)
# git clone https://github.com/masyagin1998/bio-alignment
# make
# --algo hb # hirschberg
# --algo nw # needelman
# --in query.fasta,ref.fasta

# sequence_align
# pip install sequence_align
# run script

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Print header
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

{
    printf "Program\tref\tquery\talgorithm\tflags";
    printf "\telapsedTime\tuserTime\tsystemTime";
    printf "\tmaxResidentMemory\tCPU\n";
} > "$statsFileStr";

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Run tests
#  o sec-03 sub-01:
#    - Check which reference and query I am using
#  o sec-03 sub-02:
#    - alnSeq run two-bit needleman/waterman tests
#  o sec-03 sub-03:
#    - alnSeq run mem-water tests without scaning (slow)
#  o sec-03 sub-04:
#    - alnSeq mem-water tests with scanning (realy slow)
#  o sec-03 sub-05:
#    - alnSeq run hirschberg
#  o sec-03 sub-06:
#    - bio-align run hirschberg
#  o sec-03 sub-07:
#    - run ssw_test (complete smith waterman library)
#  o sec-03 sub-08:
#    - run blank python test & sequence_align hirschberg
#  o sec-03 sub-09:
#    - Determine if alignment is to large for some tests
#  o sec-03 sub-10:
#    - Run sequence_align and biopython (aligner) needle
#  o sec-03 sub-11:
#    - Run byte versions of alnSeq needle
#  o sec-03 sub-12:
#    - Run byte versions of alnSeq water
#  o sec-03 sub-13:
#    - Run byte versions of alnSeq water query-ref-scan
#  o sec-03 sub-14:
#    - Run bioaligment needle & emboss (needle/water)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - Check which reference and query I am using
#**********************************************************

for strRef in "$smallTestRefStr" "$midTestRefStr" "$largeTestRefStr" "$ramKillerTestRefStr"; do
# For all references
  for strQuery in "$smallTestQueryStr" "$midTestQueryStr" "$largeTestQueryStr" "$ramKillerTestQueryStr"; do
  # For all queries
    iCnt=0;

    while [[ "$iCnt" -lt "$numRepI" ]]; do
    # While I have replicates to run
    
       iCnt="$((iCnt + 1))";

       if [[ "$strRef" == "$smallTestRefStr" ]]; then
         refInStr="small";
       elif [[ "$strRef" == "$midTestRefStr" ]]; then
         refInStr="mid";
       elif [[ "$strRef" == "$largeTestRefStr" ]]; then
         refInStr="large";
       elif [[ "$strRef" == "$ramKillerTestRefStr" ]]; then
         refInStr="ramKiller";
       else
         continue; # null entry??
       fi # done

       if [[ "$strQuery" == "$smallTestQueryStr" ]]; then
         queryInStr="small";
       elif [[ "$strQuery" == "$midTestQueryStr" ]]; then
         queryInStr="mid";
       elif [[ "$strQuery" == "$largeTestQueryStr" ]]; then
         queryInStr="large";
       elif [[ "$strQuery" == "$ramKillerTestQueryStr" ]]; then
         queryInStr="ramKiller";
       else
         continue; # null entry??
       fi # done


     #*****************************************************
     # Sec-03 Sub-02:
     #  - alnSeq run two-bit needleman/waterman tests
     #*****************************************************

       # Two bit version of alnSeq
       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\tneedle\ttwo-bit\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqTwoBit" \
              -use-needle \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\twater\ttwo-bit\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqTwoBit" \
              -use-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tscan-water\ttwo-bit\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqTwoBit" \
       #       -query-ref-scan \
       #       -use-water \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";


     #*****************************************************
     # Sec-03 Sub-03:
     #  - alnSeq run mem-water tests without scaning (slow)
     #*****************************************************

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tmem-water\ttwo-bit\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqTwoBit" \
       #       -use-mem-water \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\tmem-water\tbyte\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqByte" \
              -use-mem-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\tmem-water\tmid\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqMid" \
              -use-mem-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\tmem-water\tfast\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqFast" \
              -use-mem-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

     #*****************************************************
     # Sec-03 Sub-04:
     #  - alnSeq mem-water tests with scanning (realy slow)
     #*****************************************************

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tscan-mem-water\ttwo-bit\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqTwoBit" \
       #       -query-ref-scan \
       #       -use-mem-water \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tscan-mem-water\tbyte\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqByte" \
       #       -query-ref-scan \
       #       -use-mem-water \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tscan-mem-water\tmid\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqMid" \
       #       -query-ref-scan \
       #       -use-mem-water \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tscan-mem-water\tfast\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqFast" \
       #       -query-ref-scan \
       #       -use-mem-water \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";

     #*****************************************************
     # Sec-03 Sub-05:
     #  - alnSeq run hirschberg
     #*****************************************************

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\thirschberg\ttwo-bit\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqTwoBit" \
       #       -use-hirschberg \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";


       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\thirschberg\tbyte\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqByte" \
              -use-hirschberg \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";


       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\thirschberg\tmid\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqMid" \
              -use-hirschberg \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\thirschberg\tfast\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqFast" \
              -use-hirschberg \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

     #*****************************************************
     # Sec-03 Sub-06:
     #  - bio-align run hirschberg
     #*****************************************************

       # bio alignment with heirschenberg
       /usr/bin/time \
           -f "bio-align\t$refInStr\t$queryInStr\thirschberg\tNA\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../bio-alignment/bin/bio-alignment" \
              --in "$strQuery,$strRef" \
              --algo hb \
            > tmp.aln;

     #*****************************************************
     # Sec-03 Sub-07:
     #  - run ssw_test (complete smith waterman library)
     #*****************************************************

       /usr/bin/time \
           -f "ssw_test\t$refInStr\t$queryInStr\twater\tNA\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../Complete-Striped-Smith-Waterman-Library/src/ssw_test" \
             -c \
             "$strQuery" \
             "$strRef" \
           > "tmp.aln";

     #*****************************************************
     # Sec-03 Sub-08:
     #  - run blank python test & sequence_align hirschberg
     #*****************************************************

       # Python libraries
       # Load up the cache for python
       python3 "$scriptDirStr/../scripts/alignment-python-blankLoad.py" \
          "$strRef" \
          "$strQuery";
       python3 "$scriptDirStr/../scripts/alignment-python-blankLoad.py" \
          "$strRef" \
          "$strQuery";
       python3 "$scriptDirStr/../scripts/alignment-python-blankLoad.py" \
          "$strRef" \
          "$strQuery";

       # Get the no alignment times
       # load in sequence / print out / memory of loadin
       /usr/bin/time \
          -f "python\t$refInStr\t$queryInStr\tblank\tNA\t%e\t%U\t%S\t%M\t%P" \
          -o "$statsFileStr" \
          -a \
          python3 "$scriptDirStr/../scripts/alignment-python-blankLoad.py" \
          "$strRef" \
          "$strQuery";

       /usr/bin/time \
          -f "sequence_align\t$refInStr\t$queryInStr\thirschberg\tNA\t%e\t%U\t%S\t%M\t%P" \
          -o "$statsFileStr" \
          -a \
          python3 "$scriptDirStr/../scripts/sequence_align-hirsch.py" \
          "$strRef" \
          "$strQuery";

     #*****************************************************
     # Sec-03 Sub-09:
     #  - Determine if alignment is to large for some tests
     #*****************************************************

       # The byte matrix can not support the ram killer alignment
       if [[ "$strRef" == "$ramKillerTestRefStr" ]]; then
         if [[ "$strQuery" != "$smallTestQueryStr" ]]; then
           continue;
         fi # if alignment is to large
       fi # If the other probrams would run out of memory

       if [[ "$strQuery" == "$ramKillerTestQueryStr" ]]; then
         if [[ "$strRef" != "$smallTestRefStr" ]]; then
           continue;
         fi # if alignment is to large
       fi # If the other probrams would run out of memory

     #*****************************************************
     # Sec-03 Sub-10:
     #  - Run sequence_align and biopython (aligner) needle
     #*****************************************************

       /usr/bin/time \
          -f "sequence_align\t$refInStr\t$queryInStr\tneedle\tNA\t%e\t%U\t%S\t%M\t%P" \
          -o "$statsFileStr" \
          -a \
          python3 "$scriptDirStr/../scripts/sequence_align-needle.py" \
          "$strRef" \
          "$strQuery";

       /usr/bin/time \
          -f "biopython\t$refInStr\t$queryInStr\tneedle\tNA\t%e\t%U\t%S\t%M\t%P" \
          -o "$statsFileStr" \
          -a \
          python3 "$scriptDirStr/../scripts/biopython-pairwiseAligner.py" \
          "$strRef" \
          "$strQuery";

     #*****************************************************
     # Sec-03 Sub-11:
     #  - Run byte versions of alnSeq needle
     #*****************************************************

       # Byte version of alnSeq
       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\tneedle\tbyte\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqByte" \
              -use-needle \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       # Less options, byte matrix version
       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\tneedle\tmid\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqMid" \
              -use-needle \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       # no options, byte matrix and fastest
       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\tneedle\tfast\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqFast" \
              -use-needle \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

     #*****************************************************
     # Sec-03 Sub-12:
     #  - Run byte versions of alnSeq water
     #*****************************************************

       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\twater\tbyte\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqByte" \
              -use-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\twater\tmid\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqMid" \
              -use-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeq\t$refInStr\t$queryInStr\twater\tfast\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqFast" \
              -use-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

     #*****************************************************
     # Sec-03 Sub-13:
     #  - Run byte versions of alnSeq water query-ref-scan
     #*****************************************************

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tscan\tbyte\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqByte" \
       #       -query-ref-scan \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tscan\tmid\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqMid" \
       #       -query-ref-scan \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";

       #/usr/bin/time \
       #    -f "alnSeq\t$refInStr\t$queryInStr\tscan\tfast\t%e\t%U\t%S\t%M\t%P" \
       #    -o "$statsFileStr" \
       #    -a \
       #    "$scriptDirStr/../alnSeqFast" \
       #       -query-ref-scan \
       #       -query "$strQuery" \
       #       -ref "$strRef" \
       #       -out "tmp.aln";
    
     #*****************************************************
     # Sec-03 Sub-14:
     #  - Run bioaligment needle & emboss (needle/water)
     #*****************************************************

       # bio-alignment Regular
       /usr/bin/time \
           -f "bio-align\t$refInStr\t$queryInStr\tneedle\tNA\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../bio-alignment/bin/bio-alignment" \
              --in "$strQuery,$strRef" \
              --algo nw \
              --out tmp.aln;
    
       # emboss's needle
       /usr/bin/time \
           -f "emboss\t$refInStr\t$queryInStr\tneedle\tNA\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           needle \
              -bsequence "$strQuery" \
              -asequence "$strRef" \
          -gapextend -4 \
          -gapopen -1 \
              -out "tmp.aln";
    
       /usr/bin/time \
           -f "emboss\t$refInStr\t$queryInStr\twater\tNA\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           water \
              -bsequence "$strQuery" \
              -asequence "$strRef" \
          -gapextend -4 \
          -gapopen -1 \
              -out "tmp.aln";

    done # While I have replicates to run
  done  # For all queries
done # For all references
