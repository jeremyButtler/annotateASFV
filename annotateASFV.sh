#!/usr/bin/bash

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# anotateASFV.sh SOF: Start Of File
#  - Annotates an ASFV genome
#  o sec-01:
#    - Variable declerations
#  o sec-02:
#    - Check user input
#  o sec-03:
#    - Annotation
#  o sec-04:
#    - Remove temporary files, translate sequences and
#      detect if reading frames are complete
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###########################################################
# Name: anotateASFV.sh
# Call: bash anotateASFV.sh genes genes.fa genome.fa prefix
# Use:
#  - Annotates an ASFV genome
#  - This should be used as a guide, not a solution. Make
#    sure you check all output. There will be some mistakes
# Input:
#  - genes.fasta (\$1): [Required]
#    o Fasta file with genes to annotate with. Each gene
#      should start with atg and end with a stop codon.
#      - For non-coding regions ignore the error messages
#  - genome.fasta (\$2): [Required]
#    o Genome to annotate
#  - prefix (\$3): [ASFV]
#    o Prefix to add to output file names 
# Ouput:
#  - prefix-seq.fa
#    o Fasta file with the acid sequences. All sequences
#      start with a stop codon (so are reversed).
#      - Header: >sequenceName-start-end-forward/reverse
#      - forward = 1; reverse = 0
#  - prefix-aa.fa
#    o Fasta file with the translated amino acid sequences
#      for each annotation (ignore non-gene annotations).
#  - prefix-feat.tbl:
#    o Genbank feature table. Positions that you need to
#      select an alignment for are labled with ALT-FIX
#  - prefix-scores.tsv:
#    o tsv file with the selected alignments and possible
#      alternative alignments.
#  - prefix-map.aln:
#    o File with the individual alignments for each feature
#  - prefix-error.tsv:
#    o tsv file with all the features that need to be fixed
# Requires:
#  - EMBOSS
#  - standard unix programs. You already have these unless
#    your using windows. In which case use the linux shell.
###########################################################

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Variable declerations
#  o sec-01 sub-01:
#    - Variables maninpulated by the user
#  o sec-01 sub-02:
#    - Script variables
#  o sec-01 sub-03:
#    - Help message
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-01 Sub-01:
#  - Variables maninpulated by the user
#**********************************************************

dbStr="$1";  # fasta data base of genes to use
             # Each sequence is on one line
seqStr="$2"; # Fasta file with genome to annotate
prefixStr="$3"; # Prefix to call everything by

#**********************************************************
# Sec-01 Sub-02:
#  - Script variables
#**********************************************************

featHeadStr=""; # Name for feature table head (genome head)
iCnt=1;       # Counter for my loop
numGenesI=0;  # number of genes in dbStr
geneLenI=0;   # Length of gene annotating with
geneIdStr=""; # Name of gene annotating with
percI=10;     # length - length/percI is used to find the
              # min score for alternative alignments
scriptDirStr="$(dirname "$0")";

#**********************************************************
# Sec-01 Sub-03:
#  - Help message
#**********************************************************

helpStr="bash $(basename "$0") genes.fa genome.fa prefix
 Use:
  - Annotates an ASFV genome
  - This should be used as a guide, not a solution. Make
    sure you check all output. There will be some mistakes.
 Input:
  - genes.fasta (\$1): [Required]
    o Fasta file with genes to annotate with. Each gene
      should start with atg and end with a stop codon.
      - For non-coding regions ignore the error messages
  - genome.fasta (\$2): [Required]
    o Genome to annotate
  - prefix (\$3): [ASFV]
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
  - EMBOSS
  - standard unix programs. You already have these unless
    your using windows. In which case use the linux shell.
"; # Help message

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

# Check if the user requested the help message

if [[ "$1" == "-h" || "$1" == "--h" || $1 == "help" ]];then
   printf "%s\n" "$helpStr";
   exit;
fi

if [[ "$2" == "-h" || "$2" == "--h" || $2 == "help" ]];then
   printf "%s\n" "$helpStr";
   exit;
fi

if [[ "$3" == "-h" || "$3" == "--h" || $3 == "help" ]];then
   printf "%s\n" "$helpStr";
   exit;
fi

if [[ "$4" == "-h" || "$4" == "--h" || $4 == "help" ]];then
   printf "%s\n" "$helpStr";
   exit;
fi

# Check user input

if [[ "$prefixStr" == "" ]]; then prefixStr="ASFV"; fi

if [[ ! -f "$seqStr" ]]; then
   printf\
      "could not open genome to annotate (%s).\n"\
      "$seqStr";
   exit;
fi # If an invalid genome was input

if [[ ! -f "$dbStr" ]]; then
  {
     printf "could not fasta file with genes to annotate";
     printf " with (%s)\n" "$seqStr";
  }

  exit;
fi # If an invalid deta base was input

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - Annotation
#  o sec-03 sub-01:
#    - Find number of genes & print out file headers
#  o sec-03 sub-02:
#    - Find the target gene length, name, & report to user
#  o sec-03 sub-03:
#    - Extract and align the forward gene sequence
#  o sec-03 sub-04:
#    - Get the alternative alignments for the forward gene
#  o sec-03 sub-05:
#    - Extract and align the reverse gene sequence
#  o sec-03 sub-06:
#    - Get the alternative alignments for the reverse gene
#  o sec-03 sub-07:
#    - Find out if the forward or reverse score was better
#  o sec-03 sub-08:
#    - If both scores were equal, print out reverse as an
#      alternative entry.
#  o sec-03 sub-09:
#    - Print out the best scoring sequence
#  o sec-03 sub-10:
#    - Print out the feature table entry (best sequence)
#  o sec-03 sub-11:
#    - Save the best alignment and move to the next gene
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - Find number of genes & print out file headers
#**********************************************************

featHeadStr="$(sed -n '/^>/{s/^>//p;q;};' "$seqStr")";
printf\
    ">Feature\t%s\n"\
    "$featHeadStr"\
  > "$prefixStr-feat.tbl";

numGenesI="$(wc -l "$dbStr" | awk '{print $1}')";

{
   printf "Feature\tDirection\tscore\tgenomeStart";
   printf "\tmappedGeneStart\tgenomeEnd\tmappedGeneEnd"
   printf "\talingedLength\tmappedGenesLength\tmaxScore\n";
} > "$prefixStr-scores.tsv";

#**********************************************************
# Sec-03 Sub-02:
#  - Find the target gene length, name, and report to user
#**********************************************************

while [[ "$iCnt" -lt "$numGenesI" ]]; do
# Loop: Map all features
   geneIdStr="$(sed -n ""$iCnt"s/^>//p;" "$dbStr")";

   geneLenI="$(\
      sed -n "$((iCnt + 1))p;" "$dbStr" |
      awk '{print length($0); exit;}'\
   )";

   printf "on %s\tlength=%i\n" "$geneIdStr" "$geneLenI";

   #*******************************************************
   # Sec-03 Sub-03:
   #  - Extract and align the forward gene sequence
   #*******************************************************

   sed \
       -n \
      "$iCnt,$((iCnt + 1))p;"\
      "$dbStr"\
    > "$prefixStr-tmpQry.fa";

   # align the forward sequence
   "$scriptDirStr/alnSeq"\
       -use-water\
       -print-aligned\
       -print-positions\
       -ref "$seqStr"\
       -query "$prefixStr-tmpQry.fa"\
     > "$prefixStr-tmpDel.aln";

   #*******************************************************
   # Sec-03 Sub-04:
   #  - Get the alternative alignments for the forward gene
   #*******************************************************

   # Get alternate positons (only consider positions with
   # a score of 90% of the maximum)
   "$scriptDirStr/alnSeqScan"\
       -query-ref-scan\
       -only-scores\
       -min-score "$((5* (geneLenI - (geneLenI/percI)) ))"\
       -ref "$seqStr"\
       -query "$prefixStr-tmpQry.fa" |
     sed -n '2,$p; # remove header' |
     sort -n -k 2 -k 5 |
     sort -r -n -k 3 |
     awk\
        -v nameStr="$geneIdStr"\
        -v lenI="$geneLenI"\
        -v maxScoreI="$((5 * geneLenI))"\
        'BEGIN{
           getline;
           scoreI=$2;
           refStartI=$3;
           qryStartI=$4;
           refEndI=$5;
           qryEndI=$6;
         };
         { # MAIN
            # Check if the other frame is in this frame
            if($3 >= refStartI && $3 <= refEndI) next;

            # this mapping does not overlap, print it out
            if(scoreI > 0)
            { # If: something was found
               printf "%s\tforward\t%i", nameStr, scoreI;
               printf "\t%i\t%i", refStartI+1, qryStartI+1;
               printf "\t%i\t%i", refEndI+1, qryEndI+1; 
               printf "\t%i", refEndI - refStartI;
               printf "\t%i\t%i\n", lenI, maxScoreI;
            } # If: something was found

            # Get the new mapping
            scoreI=$2;
            refStartI=$3;
            qryStartI=$4;
            refEndI=$5;
            qryEndI=$6; 
         }; # MAIN
         END{
            if(scoreI > 0)
            { # If: something was found
               # Print the final scores
               printf "%s\tforward\t%i", nameStr, scoreI;
               printf "\t%i\t%i", refStartI+1, qryStartI+1;
               printf "\t%i\t%i", refEndI+1, qryEndI+1; 
               printf "\t%i", refEndI - refStartI;
               printf "\t%i\t%i\n", lenI, maxScoreI;
            } # If: something was found
         };' >> "$prefixStr-scores.tsv";

   #*******************************************************
   # Sec-03 Sub-05:
   #  - Extract and align the reverse gene sequence
   #*******************************************************

   # align the reverse sequence
   revseq\
       -sequence "$prefixStr-tmpQry.fa"\
       -auto\
       -stdout\
     > "$prefixStr-tmpQryR.fa";

   "$scriptDirStr/alnSeq"\
       -use-water\
       -print-aligned\
       -print-positions\
       -ref "$seqStr"\
       -query "$prefixStr-tmpQryR.fa"\
     > "$prefixStr-tmpDel-rev.aln";

   #*******************************************************
   # Sec-03 Sub-06:
   #  - Get the alternative alignments for the reverse gene
   #*******************************************************

   # Get alternate positons (only consider positions with
   # a score of 90% of the maximum)
   "$scriptDirStr/alnSeqScan"\
       -query-ref-scan\
       -only-scores\
       -min-score "$((5* (geneLenI - (geneLenI/percI)) ))"\
       -ref "$seqStr"\
       -query "$prefixStr-tmpQryR.fa" |
     sed -n '2,$p; # remove header' |
     sort -n -k 2 -k 5 |
     sort -r -n -k 3 |
     awk\
        -v nameStr="$geneIdStr"\
        -v lenI="$geneLenI"\
        -v maxScoreI="$((5 * geneLenI))"\
        'BEGIN{
           getline;
           scoreI=$2;
           refStartI=$3;
           qryStartI=$4;
           refEndI=$5;
           qryEndI=$6;
         };
         { # MAIN
            # Check if the other frame is in this frame
            if($3 >= refStartI && $3 <= refEndI) next;

            if(scoreI > 0)
            { # If: something was found
               # mapping does not overlap, print it out
               printf "%s\tbackwards\t%i", nameStr, scoreI;
               printf "\t%i\t%i", refEndI+1, qryEndI+1;
               printf "\t%i\t%i", refStartI+1, qryStartI+1; 
               printf "\t%i", refEndI - refStartI;
               printf "\t%i\t%i\n", lenI, maxScoreI;
            } # If: something was found

            # Get the new mapping
            scoreI=$2;
            refStartI=$3;
            qryStartI=$4;
            refEndI=$5;
            qryEndI=$6; 
         }; # MAIN
         END{
            if(scoreI > 0)
            { # If: something was found
               # Print out the final scores
               printf "%s\tbackwards\t%i", nameStr, scoreI;
               printf "\t%i\t%i", refEndI+1, qryEndI+1;
               printf "\t%i\t%i", refStartI+1, qryStartI+1; 
               printf "\t%i", refEndI - refStartI;
               printf "\t%i\t%i\n", lenI, maxScoreI;
            } # If: something was found
         };' >> "$prefixStr-scores.tsv";

   rm "$prefixStr-tmpQry"*.fa;

   #*******************************************************
   # Sec-03 Sub-07:
   #  - Find out if the forward or reverse score was better
   #*******************************************************

   # Find which direction had the best score
   forScoreI="$(\
      sed -n '17s/.*: //p' < "$prefixStr-tmpDel.aln"\
   )";

   revScoreI="$(\
      sed -n '17s/.*: //p' < "$prefixStr-tmpDel-rev.aln"\
   )";

   if [[ "$forScoreI" -gt "$revScoreI" ]]; then
      alnStr="$prefixStr-tmpDel.aln";
      forBl=1;

   elif [[ "$forScoreI" -lt "$revScoreI" ]]; then
      alnStr="$prefixStr-tmpDel-rev.aln";
      forBl=0;

   #*******************************************************
   # Sec-03 Sub-08:
   #  - If both scores were equal, print out reverse as an
   #    alternative entry.
   #*******************************************************

   else
   # Else if: There are two possible locations
      alnStr="$prefixStr-tmpDel.aln";
      forBl=1;
      # I need to print out the alternate feature and seq

      # Add the alternate location to my files
      awk\
         -v dirBl="$forBl"\
         -v tblStr="$prefixStr-feat.tbl"\
         '{
            # Get the query name (gene mapping to genome)
            if(NR == 2){headStr = $4; next;}; # gene name

            # Get the reference start and end (full genome)
            if(NR == 5)
            { # If: I am on the reference mapping position
              # Get the positon of the gene on the
              # reference. This is a reverse complement,
              # so the starting locations are reversed.
              startI = $8; # The end (first base in frame)
              endI = $6;   # The start (last base in frame)
 
              printf "%i\t%i", $6, $8 >> tblStr;

              if(headStr ~ /repeat/)
              { # If: this was a repeat region
               printf "\trepeat_region\tALT-FIX\n">>tblStr;
               printf "\t\t\trpt_type\tinverted\n">>tblStr;
               printf "\t\t\trpt_type\tterminal\n">>tblStr;
              } # If: this was a repeat region

              else
              { # Else: this was a gene
               printf "\tgene\tALT-FIX\n" >> tblStr;
               printf "\t\t\tgene\t%s\n",headStr >> tblStr;
               printf "%i\t%i", startI, endI >> tblStr;
               printf "\tCDS\tALT-FIX\n" >> tblStr;
               printf "\t\t\tproduct" >> tblStr;
               printf "\t%s\n", headStr >> tblStr;
              } # Else: this was a gene
            }; # If: I am on the reference mapping position
            
            if($1 == "Ref:"){seqStr = seqStr $3; next;};
          }; # MAIN
          END{
            printf ">%s-%i", headStr, startI;
            printf "-%i-%i-ALT-FIX\n", endI, dirBl;
            printf "%s\n\n", seqStr;
          }; # Print out the sequence
         ' < "$prefixStr-tmpDel-rev.aln"\
         > "$prefixStr-tmpDel-rev.fa";

         # Set up the reverse sequence for the correct
         # reading frame
         revseq\
             -sequence "$prefixStr-tmpDel-rev.fa"\
             -auto\
             -stdout\
           >> "$prefixStr-seq.fa";

         printf "\n\n" >> "$prefixStr-seq.fa";
         rm "$prefixStr-tmpDel-rev.fa";

         cat\
             "$prefixStr-tmpDel-rev.aln"\
           >> "$prefixStr-map.aln";

         printf "\n\n" >> "$prefixStr-map.aln";
   # Else: both scores were equal
   fi; # Check which alignment is better

   #*******************************************************
   # Sec-03 Sub-09:
   #  - Print out the best scoring sequence
   #*******************************************************

   awk\
      -v dirBl="$forBl"\
      '{
         # Get the query name (gene mapping to genome)
         if(NR == 2){headStr = $4; next;}; # Header entry

         # Get the reference start and end (full genome)
         if(NR == 5){startI = $6; endI = $8; next;};
         if($1 == "Ref:"){seqStr = seqStr $3; next;};
       }; # MAIN
       END{
         printf ">%s-%i-%i-%i\n",headStr,startI,endI,dirBl;
         printf "%s\n\n", seqStr;
       }; # Print out the sequence
    ' < "$alnStr"\
    > "$prefixStr-tmpDel.fa";
 
   # Check if this is a reverse complement alignment
   # If so, then I need to reverse the sequence to get the
   # correct reading frame
   if [[ "$forBl" == 1 ]]; then
      cat "$prefixStr-tmpDel.fa" >> "$prefixStr-seq.fa";
   else
      revseq\
          -sequence "$prefixStr-tmpDel.fa"\
          -auto\
          -stdout\
         >> "$prefixStr-seq.fa";
   fi

   rm "$prefixStr-tmpDel.fa";

   #*******************************************************
   # Sec-03 Sub-10:
   #  - Print out the feature table entry (best sequence)
   #*******************************************************

   awk\
       -v dirBl="$forBl"\
       '{
           getline;
           gsub(/-/, " ", $4);
           headStr = $4;
           getline;   # Move to query position
           getline;   # Move to ref name
           getline;   # Move to ref position
           startI = $6;
           endI = $8;

           # Check if the reading frame is backwards
           if(dirBl == 1) printf "%i\t%i", startI, endI;
           else printf "%i\t%i", endI, startI;

           if(headStr ~ /repeat/)
           { # If: this is not a gene
              printf "\trepeat_region\n";
              printf "\t\t\trpt_type\tinverted\n";
              printf "\t\t\trpt_type\tterminal\n";
              exit;
           } # If: this is not a gene
           
           printf "\tgene\n";
           printf "\t\t\tgene\t%s\n", headStr;

           if(dirBl == 1) printf "%i\t%i", startI, endI;
           else printf "%i\t%i", endI, startI;

           printf "\tCDS\n\t\t\tproduct\t%s\n", headStr;

           exit;
        }; # MAIN
      ' < "$alnStr"\
    >> "$prefixStr-feat.tbl";

   #*******************************************************
   # Sec-03 Sub-11:
   #  - Save the best alignment and move to the next gene
   #*******************************************************

   #  save the aligment
   cat "$alnStr" >> "$prefixStr-map.aln";
   printf "\n\n" >> "$prefixStr-map.aln";

   iCnt="$((iCnt + 2))";
done # Loop: Map all features

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-04:
#  - Remove temporary files, translate sequences and
#    detect if reading frames are complete
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

rm "$prefixStr-tmpDel.aln" "$prefixStr-tmpDel-rev.aln";

transeq\
    -sequence "$prefixStr-seq.fa"\
    -auto\
    -stdout\
  > "$prefixStr-aa.fa";

awk\
    -f "$scriptDirStr/detectIncompleteORF.awk"\
    < "$prefixStr-aa.fa"\
 > "$prefixStr-errors.tsv";

awk '{
    if(NR == 1) next;   # Skip the header
    gsub(/-/, " ", $1);
    gsub(/-/, " ", $2);
    print $1 " (" $4 " to " $5 ") " "has " $2;
  }'\
  < "$prefixStr-errors.tsv";

exit;
