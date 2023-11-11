# git clone https://github.com/mengyao/Complete-Striped-Smith-Waterman-Library
# cd Complete-Striped-Smith-Waterman-Library/src
# make

# git clone https://github.com/masyagin1998/bio-alignment
# cd bio-alignment
# make

# Does not work on openbsd
# git clone --recursive https://github.com/noporpoise/seq-align
# cd seq-align
# make

# git clone https://github.com/OryzaBioinformatics/emboss
# cd emboss
# ./configure
# Had issue with missing files

scriptDirStr="$(dirname "$0")";
smallTestQueryStr="$scriptDirStr/../genomes/Small-query.fasta";
smallTestRefStr="$scriptDirStr/../genomes/Small-ref.fasta";

midTestQueryStr="$scriptDirStr/../genomes/Mid-query.fasta";
midTestRefStr="$scriptDirStr/../genomes/Mid-reference.fasta";

largeTestQueryStr="$scriptDirStr/../genomes/Large-query.fasta";
largeTestRefStr="$scriptDirStr/../genomes/Large-reference.fasta";

ramKillerTestQueryStr="$scriptDirStr/../genomes/RamKiller-query.fasta";
ramKillerTestRefStr="$scriptDirStr/../genomes/RamKiller-ref.fasta";

numRepI=4;
statsFileStr="$scriptDirStr/20230722-alnSeq-bench-stats.tsv";
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

{
    printf "Program\ttest\telapsedTime\tuserTime\tsystemTime";
    printf "\tmaxResidentMemory\tCPU\n";
} > "$statsFileStr";

for strRef in "$smallTestRefStr" "$midTestRefStr" "$largeTestRefStr" "$ramKillerTestRefStr"; do
# For all references
  for strQuery in "$smallTestQueryStr" "$midTestQueryStr" "$largeTestQueryStr" "$ramKillerTestQueryStr"; do
  # For all queries
    #printf "%s\t%s\n" "$strRef" "$strQuery"; continue;
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

       #printf "%s\t%s\n" "$refInStr" "$queryInStr"; continue;

       # alnigSeq oFast
       /usr/bin/time \
           -f "alnSeqOFast\t$refInStr\t$queryInStr\tneedle\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqOfast" \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeqOFast\t$refInStr\t$queryInStr\twater\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqOfast" \
              -use-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";
    
       /usr/bin/time \
           -f "alnSeqO3\t$refInStr\t$queryInStr\tneedle\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqO3" \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeqO3\t$refInStr\t$queryInStr\twater\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqO3" \
              -use-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";
    
       # alnigSeq O2
       /usr/bin/time \
           -f "alnSeqO2\t$refInStr\t$queryInStr\tneedle\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqO2" \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       # alnigSeq O2
       /usr/bin/time \
           -f "alnSeqO2\t$refInStr\t$queryInStr\twater\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqO2" \
              -use-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";
    
       /usr/bin/time \
           -f "alnSeqO0\t$refInStr\t$queryInStr\tneedle\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqO0" \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       /usr/bin/time \
           -f "alnSeqO0\t$refInStr\t$queryInStr\twater\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../alnSeqO0" \
              -use-water \
              -query "$strQuery" \
              -ref "$strRef" \
              -out "tmp.aln";

       # bio alignment with heirschenberg
       /usr/bin/time \
           -f "bio-align\t$refInStr\t$queryInStr\thb\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../bio-alignment/bin/bio-alignment" \
              --in "$strQuery,$smallTestRefStr" \
    	  --algo hb \
              > tmp.aln;

       /usr/bin/time \
           -f "ssw_test\t$refInStr\t$queryInStr\twater\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../Complete-Striped-Smith-Waterman-Library/src/ssw_test" -c \
             "$strQuery" \
             "$strRef" \
           > "tmp.aln";
    
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

       # bio-alignment Regular
       /usr/bin/time \
           -f "bio-align\t$refInStr\t$queryInStr\tneedle\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../bio-alignment/bin/bio-alignment" \
              --in "$strQuery,$strRef" \
              --algo nw \
              --out tmp.aln;
    
       /usr/bin/time \
           -f "bio-align\t$refInStr\t$queryInStr\twater\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../bio-alignment/bin/bio-alignment" \
              --in "$strQuery,$strRef" \
    	  --algo sw \
          --out tmp.aln;
    
       # seq-align
       /usr/bin/time \
           -f "seq-align\t$refInStr\t$queryInStr\tneedle\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../seq-align/bin/needleman_wunsch" --files "$strQuery" "$strRef" > tmp.aln;
    
       /usr/bin/time \
           -f "seq-align\t$refInStr\t$queryInStr\twater\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           "$scriptDirStr/../seq-align/bin/smith_waterman" --maxhits 1 --files "$strQuery" "$strRef" > tmp.aln;
    
       # emboss's needle
       /usr/bin/time \
           -f "emboss\t$refInStr\t$queryInStr\tneedle\t%e\t%U\t%S\t%M\t%P" \
           -o "$statsFileStr" \
           -a \
           needle \
              -bsequence "$strQuery" \
              -asequence "$strRef" \
    	  -gapextend -4 \
    	  -gapopen -1 \
              -out "tmp.aln";
    
       /usr/bin/time \
           -f "emboss\t$refInStr\t$queryInStr\twater\t%e\t%U\t%S\t%M\t%P" \
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

rm "tmp.aln";
exit;
