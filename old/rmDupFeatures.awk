{ # MAIN
   if($1 ~ /^>/)
   { # If: we are starting a new sequence
      if(geneNameStr == "")
      { # If: this is the first sequence
         sub(/.*="*/, "", $0);
         sub(/"$/, "", $0);
         gsub(/[ _]/, "-", $0); # I do not like spaces or _
         geneNameStr = $0;
         next;
      } # If: this is the first sequence

      # Remove everything but the feature/gene name
      sub(/.*="*/, "", $0);
      sub(/"$/, "", $0);
      gsub(/[ _]/, "-", $0); # I do not like spaces or _'s

      # Check if the next sequence is a duplicate
      if($0 == geneNameStr) {dupBl = 1; next;};

      # Print out the old sequence
      printf ">%s\n%s\n", geneNameStr, seqStr;

      # Prepare to grab the next sequence
      geneNameStr = $0;
      dupBl = 0;
      seqStr = "";
      next;
   } # If: we are starting a new sequence

   if(dupBl == 1){next;};

   seqStr = seqStr $0;
} # MAIN

# print out the final sequence
END{printf ">%s\n%s\n", geneNameStr, seqStr;};
