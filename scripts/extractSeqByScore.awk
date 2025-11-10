#!/usr/bin/awk -f

##########################################################
# extractSeqByScore.awk
#   - extracts gene sequences from a prefix-highScores.tsv
#     or prefix-scores.tsv file output by annotateASFV.sh
# Input:
#   genome.fasta:
#     o reference fasta file used with annotateASFV
#     o would be genome you are annotating
#   <file>/stdin:
#     o prefix-highScores.tsv or prefix-scores.tsv file
#       from annnotateASFV
# Output:
#   - stdout:
#     o gene/feature sequences in a fasta file
##########################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# extractSeqByScore SOF: Start Of File
#   o sec01:
#     - get reference name
#   o sec02:
#     - get reference sequence (in one line)
#   o sec03:
#     - extract and print extracted sequence
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> # Sec01:
# Sec01:
#   - get reference name
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

BEGIN{
   getline;
   fileStr = FILENAME;
   refIdStr = $0;
   sub(/^>/, "", refIdStr);
   gsub(/[ \t].*/, "", refIdStr);
};

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> # Sec01:
# Sec02:
#   - get reference sequence (in one line)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

{ # MAIN
   if(fileStr == FILENAME)
      wholeSeqStr = wholeSeqStr $0;

   #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> # Sec01:
   # Sec03:
   #   - extract and print extracted sequence
   #   o sec03 sub01:
   #     - extract sequence and print first part of header
   #   o sec03 sub02:
   #     - add direction to header and if needed reverse
   #       complement the extracted sequence
   #   o sec03 sub03:
   #     - print rest of header and teh sequence
   #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

   #******************************************************
   # Sec03 Sub01:
   #   - extract sequence and print first part of header
   #******************************************************

   else
   { # Else: on the file entry
      if($1 == "gene_id")
         next; # skip the header

      seqStr = substr(wholeSeqStr, $5, $10);
      printf ">%s-%s-%s-%s", refIdStr, $1, $5, $6;

      #***************************************************
      # Sec03 Sub02:
      #   - add direction to header and if needed reverse
      #     complement the extracted sequence
      #***************************************************

      if($2 == "forward")
         printf "-F";
      else
      { # Else: need to revese complement
         printf "-R";

         gsub(//, " ", seqStr);
         lenSI = split(seqStr, seqArySC, " ");
         seqStr = "";

         for(siNt = lenSI; siNt > 0; --siNt)
         { # Loop: reverse complement
            if(seqArySC[siNt] == "A")
               seqStr = seqStr "T"
            else if(seqArySC[siNt] == "a")
               seqStr = seqStr "t"
            else if(seqArySC[siNt] == "T")
               seqStr = seqStr "A"
            else if(seqArySC[siNt] == "t")
               seqStr = seqStr "a"
            else if(seqArySC[siNt] == "G")
               seqStr = seqStr "C"
            else if(seqArySC[siNt] == "g")
               seqStr = seqStr "c"
            else if(seqArySC[siNt] == "C")
               seqStr = seqStr "G"
            else if(seqArySC[siNt] == "c")
               seqStr = seqStr "g"

            else if(seqArySC[siNt] == "M")
               seqStr = seqStr "K"
            else if(seqArySC[siNt] == "m")
               seqStr = seqStr "k"
            else if(seqArySC[siNt] == "K")
               seqStr = seqStr "M"
            else if(seqArySC[siNt] == "k")
               seqStr = seqStr "m"
            else if(seqArySC[siNt] == "R")
               seqStr = seqStr "Y"
            else if(seqArySC[siNt] == "r")
               seqStr = seqStr "y"
            else if(seqArySC[siNt] == "Y")
               seqStr = seqStr "R"
            else if(seqArySC[siNt] == "y")
               seqStr = seqStr "r"

            else if(seqArySC[siNt] == "B")
               seqStr = seqStr "V"
            else if(seqArySC[siNt] == "b")
               seqStr = seqStr "v"
            else if(seqArySC[siNt] == "V")
               seqStr = seqStr "B"
            else if(seqArySC[siNt] == "v")
               seqStr = seqStr "b"
            else if(seqArySC[siNt] == "H")
               seqStr = seqStr "D"
            else if(seqArySC[siNt] == "h")
               seqStr = seqStr "d"
            else if(seqArySC[siNt] == "D")
               seqStr = seqStr "h"
            else if(seqArySC[siNt] == "D")
               seqStr = seqStr "h"
         } # Loop: reverse complement
      } # Else: need to revese complement

      #***************************************************
      # Sec03 Sub03:
      #   - print rest of header and teh sequence
      #***************************************************

      if($1 ~ /non-coding/)
      { # If: non-coding entry*/
         sub(/non-coding/, "", $1);
         printf "\tnon-coding";
      } # If: non-coding entry*/

      else
         printf "\tcoding";

      printf "\tid=%s", $1;
      printf "\tstart=%s", $5;
      printf "\tend=%s", $6;
      printf "\tscore=%s", $3;
      printf "\tdirection=%s", $2;
      printf "\tlength=%s\n", $10;
      printf "%s\n", seqStr;
   } # Else: on the file entry
} # MAIN
