#!/usr/bin/sh

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# mapGene.sh:
#   - maps a gene to a reference sequence or a consensus
#   o sec01:
#     - variable declarations
#   o sec02:
#     - get and check user input
#   o sec03:
#     - print header if no older file
#   o sec04:
#     - map gene to reference and filter hits
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec01:
#   - variable declarations
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

refStr="";
geneStr="";
dirStr="forward";
prefixStr="out";
fastBl=1;
scriptDirStr="$(dirname "$0")";

helpStr="$(basename "$0") -ref reference.fa -gene gene.fa
   - uses memwaterScan to get alternative mappings, then
     filters out the duplicate mappings
Input:
   -ref reference.fa: fasta file with the reference or
      consensus sequence to get mappings for
   -gene gene.fa: fasta file with gene to map to the
      reference or consensus
   -rev: in output scoring file, use reverse instead
      of forwards
   -prefix $prefixStr: prefix to call the output file
   -fast yes: use the faster kmer scaning method to hunt
      for genes (disable with -slow)
Output:
   - prefix-scores.tsv: has filtered mappings for the
     gene to the reference
"

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec02:
#   - get and check user input
#   o sec02 sub01:
#     - get user input
#   o sec02 sub02:
#     - check user input
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#*********************************************************
# Sec02 Sub01:
#   - get user input
#*********************************************************

while [ $# -gt 0 ];
do   # Loop: get user input
   if [ "$1" = "-ref" ]; then
      shift;
      refStr="$1";
   elif [ "$1" = "-gene" ]; then
      shift;
      geneStr="$1";

   elif [ "$1" = "-prefix" ]; then
      shift;
      prefixStr="$1";
   elif [ "$1" = "-out" ]; then
      shift;
      prefixStr="$1";
   elif [ "$1" = "-o" ]; then
      shift;
      prefixStr="$1";
   elif [ "$1" = "-p" ]; then
      shift;
      prefixStr="$1";
   elif [ "$1" = "-fast" ]; then
      fastBl=1;
   elif [ "$1" = "-slow" ]; then
      fastBl=0;

   elif [ "$1" = "-r" ]; then
      dirStr="reverse";
   elif [ "$1" = "-rev" ]; then
      dirStr="reverse";
   elif [ "$1" = "-reverse" ]; then
      dirStr="reverse";

   elif [ "$1" = "-f" ]; then
      dirStr="forward";
   elif [ "$1" = "-for" ]; then
      dirStr="forward";
   elif [ "$1" = "-forward" ]; then
      dirStr="forward";

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
      printf "%s is not recognized\n" "$1"; exit;
   fi;

   shift; # move to the next argument
done # Loop: get user input

#*********************************************************
# Sec02 Sub02:
#   - check user input
#*********************************************************

if [ ! -f "$geneStr" ]; then
   printf "could not open -gene %s\n" "$geneStr"; exit;
fi;

if [ ! -f "$refStr" ]; then
   printf "could not open -ref %s\n" "$refStr"; exit;
fi;

if [ "$prefixStr" = "" ]; then
   prefixStr="out";
fi;

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec03:
#   - print header if no older file
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [ ! -f "$prefixStr-scores.tsv" ]; then
   { # print header to file
      printf "gene_id\tgene_direction\tscore\tmax_score";
      printf "\tref_start\tref_end";
      printf "\tqry_start\tqry_end\tqry_len\taln_len\n";
   } > "$prefixStr-scores.tsv"; # print header to file
fi;

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec04:
#   - map gene to reference and filter hits
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

if [ "$fastBl" -lt 1 ];
then # If: using slower waterman
   "$scriptDirStr/memwaterScan" \
       -min-score 90 \
       -ref "$refStr" \
       -qry "$geneStr" |
     tail -n+2 |
     sort -n -k 3 -k 7 |
     sort -r -n -k 6 |
     awk \
        -v dirStr="$dirStr" \
        'BEGIN{
           # get the first mapping
           getline;
           nameStr=$1;
           scoreSI=$3;
           maxScoreSI=$5;
           maxLenSI=$11;
           refStartSI=$6;
           refEndSI=$7;
           qryStartSI=$9;
           qryEndSI=$10;
         };
         { # MAIN
            # check if the other frame is in this frame
            if($6 >= refStartSI && $6 <= refEndSI)
               next;
     
            # this mapping does not overlap, print it out
            if(scoreSI > 0)
            { # If: something was found
               printf "%s\t%s", nameStr, dirStr;
               printf "\t%i\t%i", scoreSI, maxScoreSI;
               printf "\t%i\t%i", refStartSI, refEndSI;
               printf "\t%i\t%i", qryStartSI, qryEndSI; 
               printf "\t%i", maxLenSI;
               printf "\t%i\n", 1 + refEndSI - refStartSI;
            } # If: something was found
     
            # Get the new mapping
            nameStr=$1;
            scoreSI=$3;
            maxScoreSI=$5;
            maxLenSI=$11;
            refStartSI=$6;
            refEndSI=$7;
            qryStartSI=$9;
            qryEndSI=$10;
         }; # MAIN
     
         END{
            if(scoreSI > 0)
            { # If: something was found
               printf "%s\t%s", nameStr, dirStr;
               printf "\t%i\t%i", scoreSI, maxScoreSI;
               printf "\t%i\t%i", refStartSI, refEndSI;
               printf "\t%i\t%i", qryStartSI, qryEndSI; 
               printf "\t%i", maxLenSI;
               printf "\t%i\n", 1 + refEndSI - refStartSI;
            } # If: something was found
         };' \
     >> "$prefixStr-scores.tsv";
else
# Else: using faster kmer scan
   "$scriptDirStr/demux" \
       -perc-score 0.90 \
       -gene "$geneStr" \
       -kmer-len 7 \
       -kmer-perc 0.60 \
       "$refStr" |
     awk \
        'BEGIN{
           getline; # get past header
           getline; # get the first mapping

           nameStr = $2;
           dirStr = $3;
           scoreSI = $4;
           maxScoreSI = $5;
           lenSI = $6;
           maxLenSI = $7;
           refStartSI = $8;
           refEndSI = $9;
           qryStartSI = "NA"; # demux does not give this
           qryEndSI = "NA";   # demux does not give this

           if(dirStr == "R")
              dirStr = "reverse";
           else
              dirStr = "forward";

           if(scoreSI < 0)
              scoreSI *= -1;
         };

         { # MAIN
            # check if the other frame is in this frame
            if($6 >= refStartSI && $6 <= refEndSI)
               next;
     
            # this mapping does not overlap, print it out
            if(scoreSI > 0)
            { # If: something was found
               printf "%s\t%s", nameStr, dirStr;
               printf "\t%i\t%i", scoreSI, maxScoreSI;
               printf "\t%i\t%i", refStartSI, refEndSI;
               printf "\t%s\t%s", qryStartSI, qryEndSI; 
               printf "\t%i\t%i\n", maxLenSI, lenSI;
            } # If: something was found
     
            nameStr  =$2;
            dirStr = $3;
            scoreSI = $4;
            maxScoreSI = $5;
            lenSI = $6;
            maxLenSI = $7;
            refStartSI = $8;
            refEndSI = $9;
            qryStartSI = "NA"; # demux does not give this
            qryEndSI = "NA";   # demux does not give this

            dirStr = $3;
            if(dirStr == "R")
               dirStr = "reverse";
            else
               dirStr = "forward";

            if(scoreSI < 0)
               scoreSI *= -1;
         }; # MAIN
     
         END{
            if(scoreSI > 0)
            { # If: something was found
               printf "%s\t%s", nameStr, dirStr;
               printf "\t%i\t%i", scoreSI, maxScoreSI;
               printf "\t%i\t%i", refStartSI, refEndSI;
               printf "\t%s\t%s", qryStartSI, qryEndSI; 
               printf "\t%i\t%i\n", maxLenSI, lenSI;
            } # If: something was found
         };' \
     >> "$prefixStr-scores.tsv";
# Else: using faster kmer scan
fi;
