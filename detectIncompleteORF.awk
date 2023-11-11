###########################################################
# Name: detectIncompleteORF
# Call: awk -f detectIncompleteORF.awk aa-seq.fasta
# Use:
#  - Detects if reading frames have a start codon and one
#    stop codon at the end (repeat is ignored)
# Input:
#   ammino-acid.fasta:
#     - File with amino acid sequences to check.
#     - The names #       in the file should end with
#       -start-end-0 or -start-end-1.
#       - -start is the first base in the alignment
#       - -end is the last base (largest nubmer) in the
#         alignment
#       - -1 is for forwards and -0 is for backwards ORFs
# Output:
#   -stdout:
#    o Prints tsv file to stdout
#    o "gene-name\tproblem\tgene-direction\tgene-start
#       \tgene-end\tproblem-base
#       \tlast-base-in-problem-codon"
###########################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# detectIncompleteORF SOF: Start Of File
#  - Detect if open reading frames were translated
#    correctly
#  o sec-01:
#     - Begin block; print out the header
#  o sec-02:
#     - Main block find problematic sequences
#  o sec-03:
#     - End block; Check and print out final error
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-01:
#  - Begin block; print out the header
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

BEGIN{
   printf "Gene\tError\tgeneDirection";
   printf "\tGeneStart\tGeneEnd\tErrorBase\tErrorEnd\n";
}

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-02:
#  - Main block find problematic sequences
#  o sec-02 sub-01:
#    - if on a header find the first stop codons position
#  o sec-02 sub-02:
#    - Check if have a complete reading frame
#  o sec-02 sub-03:
#    - Get and remove meta data (start, end, direction)
#      from the header
#  o sec-02 sub-04:
#    - Not on header, merge the sequence entry
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-02 Sub-01:
#  - if on a header find the first stop codons position
#**********************************************************

{ # MAIN
   if($1 ~ /^>/)
   { # If: this is a header
      if(seqStr != "" && headStr !~ /repeat/)
      { # If: there was a previous sequence entry
         lenSeqI = length(seqStr);
         stopIndexI = index(seqStr, "*");
         stopPosI = stopIndexI; # So I can check stopIndex
         stopPosI *= 3; # Convert to a base cordinate

         if(forBl == 1)
         { # If: this is a foward direction
            stopPosI += seqStartI;
            endPosI = stopPosI + 2;
         } # If: this is a foward direction

         else
         { # Else: this was a reverse direction
            stopPosI = seqStartI - stopPosI;
            endPosI = stopPosI - 2;
         } # Else: this was a reverse direction

         #*************************************************
         # Sec-02 Sub-02:
         #  - Check if have a complete reading frame
         #*************************************************

         if(seqStr !~ /^M/)
         { # If: there is no start codon
            printf "%s\tno-start-codon", headStr;
            printf "\t%s\t%i\t%i",dirStr,seqStartI,seqEndI;
            printf "\tNA\tNA\n";
         } # If: there is no start codon

         if(stopIndexI == 0)
         { # If: there is no stop codon
            printf "%s\tno-stop-codon", headStr;
            printf "\t%s\t%i\t%i",dirStr,seqStartI,seqEndI;
            printf "\tNA\tNA\n";
         } # If: there is no stop codon

         else if(stopIndexI < lenSeqI)
         { # Else if: there is an early stop codon
            printf "%s\tearly-stop-codon", headStr;
            printf "\t%s\t%i\t%i",dirStr,seqStartI,seqEndI;
            printf "\t%i\t%i\n", stopPosI, endPosI;
         } # Else if: there is an early stop codon
      } # If: there was a prefious sequence entry

      #****************************************************
      # Sec-02 Sub-03:
      #  - Get and remove meta data (start, end, direction)
      #    from the header
      #****************************************************

      headStr = $0;
      seqStr = "";
      sub(/-*_[0-9]$/, "", headStr); # remove reading frame

      # Remove flags I might have added in
      sub(/-*A*L*T*-*F*I*X*$/, "", headStr);
 
      # Find if this was a foward or reverse sequence
      forBl = headStr;
      sub(/.*-/, "", forBl); # Get the direction

      seqStartI = headStr;
      sub(/-[01]$/, "", seqStartI); # remove the direction

      if(forBl == 1)
      { # If: this was a forward sequence
         seqEndI = seqStartI;
         sub(/.*-/, "", seqEndI); # get end pos

         sub(/-[0-9]*$/, "", seqStartI); # remove end pos
         sub(/.*-/, "", seqStartI); # get start position
         dirStr = "Forward";
      } # If: this was a forward sequence

      # Reverse sequences start at the end
      else
      { # Else: this is a reverse sequence
         seqEndI = seqStartI;
         sub(/-[0-9]*$/, "", seqEndI); # remove start pos
         sub(/.*-/, "", seqEndI);      # get end position

         sub(/.*-/, "", seqStartI); # Start is end
         dirStr = "Reverse";
      } # Else: this is a reverse sequence

      # Remove start, end, and direction from the header
      sub(/-[0-9]*-[0-9]*-*[0-9]*$/, "", headStr);
      sub(/^>/, "", headStr);
   } # If: this is a header

   #****************************************************
   # Sec-02 Sub-04:
   #  - Not on header, merge the sequence entry
   #****************************************************

   else seqStr = seqStr $0;
} # MAIN

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec-03:
#  - End block; Check and print out final error
#  o sec-03 sub-01:
#    - Find the first stop codons position
#  o sec-03 sub-02:
#    - Check if have a complete reading frame
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

#**********************************************************
# Sec-03 Sub-01:
#  - find the first stop codons position
#**********************************************************

END{
   if(seqStr != "" && headStr !~ /repeat/)
   { # If: this is not blank or a repeat
      lenSeqI = length(seqStr);
      stopIndexI = index(seqStr, "*");
      stopPosI = stopIndexI; # So I can check stopIndex
      stopPosI *= 3; # Convert to a base cordinate

      if(forBl == 1)
      { # If: this is a foward direction
         stopPosI += seqStartI;
         endPosI = stopPosI + 2;
      } # If: this is a foward direction

      else
      { # Else: this was a reverse direction
         stopPosI = seqStartI - stopPosI;
         endPosI = stopPosI - 2;
      } # Else: this was a reverse direction

      #****************************************************
      # Sec-03 Sub-02:
      #  - Check if have a complete reading frame
      #****************************************************

      if(seqStr !~ /^M/)
      { # If: there is no start codon
         printf "%s\tno-start-codon", headStr;
         printf "\t%s\t%i\t%i",dirStr,seqStartI,seqEndI;
         printf "\tNA\tNA\n";
      } # If: there is no start codon

      if(stopIndexI == 0)
      { # If: there is no stop codon
         printf "%s\tno-stop-codon", headStr;
         printf "\t%s\t%i\t%i",dirStr,seqStartI,seqEndI;
         printf "\tNA\tNA\n";
      } # If: there is no stop codon

      else if(stopIndexI < lenSeqI)
      { # Else if: there is an early stop codon
         printf "%s\tearly-stop-codon", headStr;
         printf "\t%s\t%i\t%i",dirStr,seqStartI,seqEndI;
         printf "\t%i\t%i\n", stopPosI, endPosI;
      } # Else if: there is an early stop codon
   } # If: this is not blank or a repeat
} # END
