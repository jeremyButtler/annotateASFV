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
ramKillerTestRefStr="$scriptDirStr/../genomes/RamKiller-reference.fasta";

numRepI=10;
statsFileStr="$scriptDirStr/20231029-alnSeqPy-bench-stats.tsv"
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
#    - Run blank python test and hirschbergs
#  o sec-03 sub-03:
#    - Check if I can run the Needle/Water aligners
#  o sec-03 sub-04:
#    - Run Needle/Water alignments
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
          -f "python-alnSeq\t$refInStr\t$queryInStr\thirschberg\tNA\t%e\t%U\t%S\t%M\t%P" \
          -o "$statsFileStr" \
          -a \
       python3 "$scriptDirStr/../scripts/pyAlnSeqHirsch-test.py" \
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
     # Sec-03 Sub-03:
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
     # Sec-03 Sub-04:
     #  - Run sequence_align and biopython (aligner) needle
     #*****************************************************
       /usr/bin/time \
          -f "python-alnSeq\t$refInStr\t$queryInStr\tneedle\tNA\t%e\t%U\t%S\t%M\t%P" \
          -o "$statsFileStr" \
          -a \
       python3 "$scriptDirStr/../scripts/pyAlnSeqNeedle-test.py" \
          "$strRef" \
          "$strQuery";

     /usr/bin/time \
          -f "python-alnSeq\t$refInStr\t$queryInStr\twater\tNA\t%e\t%U\t%S\t%M\t%P" \
          -o "$statsFileStr" \
          -a \
       python3 "$scriptDirStr/../scripts/pyAlnSeqWater-test.py" \
          "$strRef" \
          "$strQuery";

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

    done # While I have replicates to run
  done  # For all queries
done # For all references
