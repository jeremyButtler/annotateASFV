#!/usr/bin/bash

###########################################################
# Name: anotateASFV.sh
#   - annotates a genome uses genes from a known sequence
# Input:
#  -gene genes.fasta: [Required]
#    o fasta file with genes to annotate with.
#      - should be full gene (start codon to end codon)
#      - for non-coding regions ignore the error messages
#  -ref genome.fasta: [Required]
#    o genome to annotate
#  -prefix out:
#    o prefix to add to output file names 
# Ouput:
#  - prefix-scores.tsv:
#    o tsv file with the scores for each genes possible
#      mappings
#  - prefix-highScores.tsv:
#    o tsv file with the selected best score(s) for each
#      gene
#  - prefix-ntSequences.fa
#    o fasta file with nucleotide sequence (coordinates
#      are from prefix-highScores)
#  - prefix-aaSequences.fa
#    o fasta file with the amino acid sequences for each
#      sequence in prefix-ntSequences.fa
#  - prefix-featureTable.tsv:
#    o Genbank feature table. Positions that you need to
#      select an alignment for are labled with
#      WARNING_ERROR
#  - prefix-featureErrors.aln:
#    o list annotations that need to be double checked
#      * dupications (make final decision)
#      * incomplete reading frames (missing start/stop)
###########################################################


#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# anotateASFV.sh SOF: Start Of File
#  - annotates a genome
#  o sec01:
#    - variable declerations
#  o sec02:
#    - get and check user input
#  o sec03:
#    - Annotation
#  o sec04:
#    - Remove temporary files, translate sequences and
#      detect if reading frames are complete
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec01:
#  - Variable declerations
#  o sec01 sub01:
#    - general variables
#  o sec01 sub02:
#    - help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec01 Sub01:
#   - general variables
#**********************************************************

genesStr="$1";  # fasta data base of genes to use
             # Each sequence is on one line
refStr="$2"; # Fasta file with genome to annotate
prefixStr="$3"; # Prefix to call everything by

iCnt=1;       # Counter for my loop
numGenesI=0;  # number of genes in genesStr
scriptDirStr="$(dirname "$0")";

#**********************************************************
# Sec01 Sub02:
#   - help message
#**********************************************************

helpStr="bash $(basename "$0") -gene genes.fa -ref genome.fa -prefix <good_name>
 Use:
  - Annotates an ASFV genome
  - This should be used as a guide, not a solution. Make
    sure you check all output. There will be some mistakes.
 Input:
  -gene genes.fasta: [Required]
    o Fasta file with genes to annotate with. Each gene
      should start with atg and end with a stop codon.
      - For non-coding regions ignore the error messages
      - use >non-coding<id> for non-coding regions
      - remove duplicate entries
  -ref genome.fasta: [Required]
    o Genome to annotate
  -prefix $prefixStr:
    o Prefix to add to output file names 
 Ouput:
  - prefix-seq.fa
    o Fasta file with the acid sequences. All sequences
      start with a stop codon (so are reversed).
      - Header: >sequenceName-start-end-forward/reverse
      - forward = 1; reverse = 0
  - prefix-aa.fa
    o Fasta file with the translated amino acid sequences
      for each annotation (ignore non-gene annotations).
  - prefix-feat.tbl:
    o Genbank feature table. Positions that you need to
      select an alignment for are labled with ALT-FIX
  - prefix-scores.tsv:
    o tsv file with the selected alignments and possible
      alternative alignments.
  - prefix-map.aln:
    o File with the individual alignments for each feature
  - prefix-error.tsv:
    o tsv file with all the features that need to be fixed
 Requires:
  - several programs with bioTools, which should be
    included in the repository
  - standard unix programs. You already have these unless
    your using windows. In which case use the linux shell.
"; # Help message

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec02:
#   - get and check user input
#   o sec02 sub01:
#     - get user input
#   o sec02 sub02:
#     - check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*********************************************************
# Sec02 Sub01:
#   - get user input
#*********************************************************

while [ $# -gt 0 ];
do   # Loop: get input
   if [ "$1" = "-ref" ]; then
      shift;
      refStr="$1";
   elif [ "$1" = "-gene" ]; then
      shift;
      genesStr="$1";
   elif [ "$1" = "-prefix" ]; then
      shift;
      prefixStr="$1";
   elif [ "$1" = "-out" ]; then
      shift;
      prefixStr="$1";
   elif [ "$1" = "-o" ]; then
      shift;
      prefixStr="$1";

   elif [ "$1" = "-h" ]; then
      printf "%s\n" "$helpStr"; exit;
   elif [ "$1" = "--h" ]; then
      printf "%s\n" "$helpStr"; exit;
   elif [ "$1" = "help" ]; then
      printf "%s\n" "$helpStr"; exit;
   elif [ "$1" = "-help" ]; then
      printf "%s\n" "$helpStr"; exit;
   elif [ "$1" = "--help" ]; then
      printf "%s\n" "$helpStr"; exit;

   else
      printf "%s is not recongnized" "$1"; exit;
   fi;

   shift;
done # Loop: get input

#*********************************************************
# Sec02 Sub02:
#   - check user input
#*********************************************************

if [ "$prefixStr" = "" ]; then prefixStr="ASFV"; fi;

if [ ! -f "$refStr" ]; then
   printf\
      "could not open genome to annotate (%s).\n"\
      "$refStr";
   exit;
fi # If an invalid genome was input

if [ ! -f "$genesStr" ]; then
  {
     printf "could not fasta file with genes to annotate";
     printf " with (%s)\n" "$refStr";
  }

  exit;
fi # If an invalid deta base was input

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec03:
#  - annotation
#  o sec03 sub01:
#    - find number genes + start loop + extract each gene
#  o sec03 sub02:
#    - get mappings for the gene
#  o sec03 sub03:
#    - get high scoring mapping(s) + merge scores + clean
#  o sec03 sub04:
#    - extract genes from genome, tanslate, and build
#      feature table
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec03 Sub01:
#  - find number genes + start loop + extract each gene
#**********************************************************

numGenesI="$(grep -c '>' "$genesStr")";
iCnt=0;
while [ "$iCnt" -lt "$numGenesI" ];
do # Loop: git hits and extract genes in reference
   iCnt="$((iCnt + 1))";

   # extract the sequence
   awk \
       -v targSI="$iCnt" \
       '
          BEGIN{seqOneSI = 0; headBl = 1;};

          { # MAIN
             if($1 ~ /^>/)
                ++seqOnSI;

             if(seqOnSI == targSI)
             { # If: On the correct sequence
                if(headBl == 1)
                { # If: on the header
                   headBl = 0;
                   sub(/>gene /, ">", $0);
                   sub(/>CDS /, ">", $0);
                   sub(/\/gene=/, "", $0);
                   gsub(/"/, "", $0);
                   gsub(/[ \t]/, "_", $0);
                   print $0;
                } # If: on the header

                else
                  print $0;
             } # If: On the correct sequence

             else if(seqOneSI > targSI)
                exit; # past the correct sequence
          } # MAIN
       ' "$genesStr" \
     > "$prefixStr-tmpQryForward.fa";

   printf "on gene %s\n" \
       "$(head -n 1 "$prefixStr-tmpQryForward.fa")";

   #*******************************************************
   # Sec03 Sub02:
   #  - get mappings for the gene
   #*******************************************************

   if [ -f "$prefixStr-scores.tsv" ]; then
      tmpPrefixStr="$prefixStr-delete-tmp";
   else
      tmpPrefixStr="$prefixStr";
   fi;

   #_____________forward_sequence_mappings________________
   sh "$scriptDirStr/mapGene.sh" \
      -ref "$refStr" \
      -gene "$prefixStr-tmpQryForward.fa" \
      -forward \
      -prefix "$tmpPrefixStr";

   #_____________reverse_sequence_mappings________________
   "$scriptDirStr/revCmp" \
       -out "$prefixStr-tmpQryReverse.fa" \
       "$prefixStr-tmpQryForward.fa";

   sh "$scriptDirStr/mapGene.sh" \
      -ref "$refStr" \
      -gene "$prefixStr-tmpQryReverse.fa" \
      -reverse \
      -prefix "$tmpPrefixStr";

   #******************************************************
   # Sec03 Sub03:
   #  - get high scoring mapping(s) + merge scores + clean
   #******************************************************

   if [ ! -f "$prefixStr-highScores.tsv" ]; then
      head -n 1 "$tmpPrefixStr-scores.tsv" \
        > "$prefixStr-highScores.tsv";
   fi;

   numLinesSI="$( \
      wc -l "$tmpPrefixStr-scores.tsv" | cut -d " " -f 1 \
   )";
   if [ "$numLinesSI" -gt 1 ];
   then # If: had matches
      # get the highscoring entries
      tail -n+2 "$tmpPrefixStr-scores.tsv" |
        sort -r -n -k 3 |
        awk '
              BEGIN{getline; print $0; highScoreSI=$3; };
              {if($3 == highScoreSI) print $0; else exit;}
            ' >> "$prefixStr-highScores.tsv";

      if [ "$tmpPrefixStr" = "$prefixStr" ]; then
         tmpPrefixStr="";
      else
         tail -n+2 "$tmpPrefixStr-scores.tsv" \
           >> "$prefixStr-scores.tsv";
         rm "$tmpPrefixStr-scores.tsv";
      fi;
   else
      if [ "$tmpPrefixStr" = "$prefixStr" ]; then
         tmpPrefixStr="";
      else
         rm "$tmpPrefixStr-scores.tsv";
      fi;
   fi; # If: had matches

   rm "$prefixStr-tmpQryReverse.fa";
   rm "$prefixStr-tmpQryForward.fa";
done # Loop: git hits and extract genes in reference

#*********************************************************
# Sec03 Sub04:
#  - extract genes from genome, tanslate, & build feature
#    table
#*********************************************************

# check if we found any genes
numLinesSI="$( \
   wc -l "$prefixStr-highScores.tsv" | cut -d " " -f 1 \
)";

if [ "$numLinesSI" -le 1 ]; then
   printf "no genes were found for -ref %s\n" "$refStr";
   exit;
fi;

# extract genes and align
awk \
    -f "$scriptDirStr/extractSeqByScore.awk" \
    "$refStr" \
    "$prefixStr-highScores.tsv" \
  > "$prefixStr-ntSequences.fa"

# translate extracted genes
"$scriptDirStr/tranSeq" \
    -orf-1 \
    -fa "$prefixStr-ntSequences.fa" \
  > "$prefixStr-aaSequences.fa"

# build feature table
awk \
    -f "$scriptDirStr/aaToFeature.awk" \
    -v prefixStr="$prefixStr" \
    "$prefixStr-aaSequences.fa" \
  > "$prefixStr-featureTable.tsv";
