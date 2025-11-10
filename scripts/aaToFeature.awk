#!/usr/bin/awk -f

##########################################################
# aaToFeature.awk
#   - converts the amino acid sequences file from
#     annotateASFV (prefix-aaSequences.fa) and converts it
#     to a Genbank feature table
# Input:
#   -v prefixStr="<prefix>":
#      o prefix to call the output error file
#   <file>/stdin:
#      o prefix-aaSequences.fa file from annotateASFV
# Output:
#   - stdout:
#     o feature table for Genbank
#   - prefix-featureErrors.tsv:
#     o tsv file with potential errors for features
##########################################################

#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
# aaToFeature.awk SOF: Start Of File
#   o sec01:
#     - setup error file header and get first aa sequence
#   o sec02:
#     - get next features and check if is a duplicate gene
#   o sec03:
#     - print out the previous feature
#   o sec04:
#     - get start, end, & coding status of next feature
#   o sec05:
#     - check if feature (next) has complete reading frame
#   o sec06:
#     - print the final feature and its errors (if any)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec01:
#   - setup error file header and get first aa sequence
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

BEGIN{
   OFS="\t";

   if(prefixStr == "")
      errFileStr = "out-featureErrors.tsv";
   else
      errFileStr = prefixStr "-featureErrors.tsv";
   printf "gene\tduplicate\tno_orf_start" > errFileStr;
   printf "\tno_orf_end\tdirection\tstart" >> errFileStr;
   printf "\tend\n" >> errFileStr;

   getline;

   # header
   sub(/-.*/, "", $1);
   sub(/^>/, "", $1);
   printf ">Feature\t%s\n", $1;

   sub("id=", "", $3);
   lastIdStr = $3;

   sub("start=", "", $4);
   startSI = $4;

   sub("end=", "", $5);
   endSI = $5;

   if($7 == "direction=reverse")
      revBl = 1;
   else
      revBl = 0;

   if($2 == "non-coding")
      nonCodingBl = 1;
   else
      nonCodingBl = 0;
}; # BEGIN

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec02:
#   - get next features id and check if is duplicate gene
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

{ # MAIN
  if($0 ~ /^>/)
  { # If: on header
    sub("id=", "", $3);
    idStr = $3;

    if(revBl)
       posStr = endSI "\t" startSI;
    else
       posStr = startSI "\t" endSI;

     if(idStr == lastIdStr)
     { # If: duplicate id

        if(warnStr != "" && ! dupIdBl)
           warnStr = warnStr "____WARNING_DUPLICATE";
        else if(! dupBl)
           warnStr = "WARNING_DUPLICATE";

        dupIdBl = 1;
        lastIdDupBl = 1;
     } # If: duplicate id

     else if(lastIdDupBl)
     { # Else If: last gene had a duplicate
        if(warnStr != "" && ! dupIdBl)
           warnStr = warnStr "____WARNING_DUPLICATE";
        else if(! dupBl)
           warnStr = "WARNING_DUPLICATE";

        lastIdDupBl = 0;
        dupIdBl = 1;
     } # Else If: last gene had a duplicate

     else
        dupIdBl = 0;

     #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     # Sec03:
     #   - print out the previous feature
     #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     # feature table entry
     if(nonCodingBl)
        printf "%s\t%s", posStr, lastIdStr;

     else 
     { # Else: is a coding (gene) feature
        printf "%s\tgene\n", posStr;
        printf "\t\t\t\tgene\t%s\n", lastIdStr;

        printf "%s\tCDS\n", posStr;
        printf "\t\t\t\tproduct\t%s\n", lastIdStr;
     } # Else: is a coding (gene) feature

     # error file warnings (and featue table warning)
     if(warnStr != "")
     { # If: printing the warning
        printf "\t\t\t\tnote\t%s\n", warnStr;

        printf "%s", lastIdStr >> errFileStr;
        if(dupIdBl)
           printf "\tTrue" >> errFileStr;
        else
           printf "\tFalse" >> errFileStr;

        if(noStartBl)
           printf "\tTrue" >> errFileStr;
        else
           printf "\tFalse" >> errFileStr;

        if(noEndBl)
           printf "\tTrue" >> errFileStr;
        else
           printf "\tFalse" >> errFileStr;
         
        if(revBl)
           printf "\treverse" >> errFileStr;
        else
           printf "\tforward" >> errFileStr;

        printf "\t%s\t%s\n", startSI, endSI >> errFileStr;
     } # If: printing the warning

     #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
     # Sec04:
     #   - get start, end, & coding status of next feature
     #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

     lastIdStr = idStr;

     sub("start=", "", $4);
     startSI = $4;

     sub("end=", "", $5);
     endSI = $5;

     if($7 == "direction=reverse")
        revBl = 1;
     else
        revBl = 0;

     if($2 == "non-coding")
        nonCodingBl = 1;
     else
        nonCodingBl = 0;

     partialORFBl = 0;
     noStartBl = 0;
     noEndBl = 0;

     if(lastIdDupBl)
        warnStr = "WARNING_DUPLICATE";
     else
        warnStr = "";
  } # If: on header

  #>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  # Sec05:
  #   - check if feature (next) has complete reading frame
  #<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  else
  { # Else: on sequence
     if($0 !~ /^[mM].*\*$/ && ! nonCodingBl)
     { # If: not a complete reading frame
        partialORFBl = 1;

        if($0 !~ /^[mM]/)
           noStartBl = 1;
        if($0 !~ /\*$/)
           noEndBl = 1;

        if(warnStr != "")
           warnStr = warnStr "____WARNING_INCOMPLETE_ORF";
        else
           warnStr = "WARNING_INCOMPLETE_ORF";
     } # If: not a complete reading frame
  } # Else: on sequence
}; # MAIN

#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
# Sec06:
#   - print the final feature and its errors (if any)
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

END{
   if(nonCodingBl)
      printf "%s\t%s", posStr, lastIdStr;

   else 
   { # Else: is a coding (gene) feature
      printf "%s\tgene\n", posStr;
      printf "\t\t\t\tgene\t%s\n", lastIdStr;

      printf "%s\tCDS\n", posStr;
      printf "\t\t\t\tproduct\t%s\n", lastIdStr;
   } # Else: is a coding (gene) feature

   if(warnStr != "")
   { # If: printing the warning
      printf "\t\t\t\tnote\t%s\n", warnStr;

      printf "%s", lastIdStr >> errFileStr;
      if(dupIdBl)
         printf "\tTrue" >> errFileStr;
      else
         printf "\tFalse" >> errFileStr;

      if(noStartBl)
         printf "\tTrue" >> errFileStr;
      else
         printf "\tFalse" >> errFileStr;

      if(noEndBl)
         printf "\tTrue" >> errFileStr;
      else
         printf "\tFalse" >> errFileStr;
       
      if(revBl)
         printf "\treverse" >> errFileStr;
      else
         printf "\tforward" >> errFileStr;

      printf "\t%s\t%s\n", startSI, endSI >> errFileStr;
   } # If: printing the warning
}; # END
