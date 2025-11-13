/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' seqStats SOF: Start Of File
'   - get stats for fx (fastq/fasta) files
'   o header:
'     - included libraries
'   o .c fun01: ln_seqStats
'     - natureal log value without using math.h (requies
'       -lm flag)
'     - TODO: fix issues with low values
'   o .c fun02: log10_seqStats
'     - log10 value without using math.h (requies -lm
'       flag)
'   o .c fun03: pow10_seqStats
'     - rases 10 to the input exponant
'   o .c fun04: qhistToMedQ_seqStats
'     - converts a q-score histogram to medain q-score
'   o .c fun05: qhistToMeanQ_seqStats
'     - converts a q-score histogram to mean q-score
'   o fun06: meanMedQ_seqStats
'     - gets the mean and median q-score for a sequence
'   o fun07: ontMeanMedQ_seqStats
'     - gets the mean and median q-score for a sequence
'       using ONT's method
'   o fun08: pReadStats_seqStats
'     - prints stats for a single read
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - included libraries
\-------------------------------------------------------*/

#ifdef PLAN9
   #include <u.h>
   #include <libc.h>
#else
   #include <stdlib.h>
#endif

#include <stdio.h>

#ifdef NO_MATH_H
#else
   #include <math.h>
#endif

#include "../genLib/endLine.h"

/*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\
! Hidden Libraries:
!   - .c  #include "../genLib/ulCp.h"
!   - .c  #include "../genLib/fileFun.h"
\%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%*/

#define def_qAdj_seqStats 33 /*convert char to q-score*/
#define def_maxQ_seqStats 94 /*highest q-score possible*/

/*-------------------------------------------------------\
| Fun01: ln_seqStats
|   - natureal log value without using math.h (requies
|     -lm flag)
'   - TODO: fix issues with low values
|   - from:
|     https://stackoverflow.com/questions/64894566/how-to-write-a-simple-logarithm-function-without-math-h
| Input:
|   - valDbl:
|     o value to find log10 for
| Output:
|   - Returns:
|     o log10 of valDbl
\-------------------------------------------------------*/
double
ln_seqStats(
   double valDbl
){
   double retDbl = 0;
   double powerDbl = 0;
   signed int roundSI = 0;
   signed int adjustSI = 0;

   #define EULER_CONST_fun03 2.718281828459045235f
      /*ln(EULER_CONST) = 1*/
   #define TAYLOR_ITERATIONS_fun03 8

   if(valDbl <= 0)
      return 0.0/0.0;

   /*move this to satisify the expansion range value
   `  of -1 < x <= 1
   */
   while(valDbl > 1.0)
   { /*Loop: remove whole logs and to taylor const min*/
      valDbl /= EULER_CONST_fun03;
      ++adjustSI;
   } /*Loop: remove whole logs and to taylor const min*/

   while(valDbl < 0.25)
   { /*Loop: move to lower talyor constant range*/
      valDbl *= EULER_CONST_fun03;
      --adjustSI;
   } /*Loop: move to lower talyor constant range*/

   /*find the decimal part of the natural log*/
   valDbl *= -1;
   powerDbl = valDbl;
   roundSI = 1;

   while(roundSI < TAYLOR_ITERATIONS_fun03)
   { /*Loop: find log2*/
      retDbl += powerDbl / roundSI;
      powerDbl *= valDbl;
      ++roundSI;

      retDbl += (powerDbl * -1) / roundSI;
      powerDbl *= valDbl;
      ++roundSI;
   }  /*Loop: find log2*/

   return retDbl + adjustSI;
} /*ln_seqStats*/

/*-------------------------------------------------------\
| Fun02: log10_seqStats
|   - log10 value without using math.h (requies -lm flag)
| Input:
|   - valDbl:
|     o value to find log10 for
| Output:
|   - Returns:
|     o log10 of valDbl
\-------------------------------------------------------*/
double
log10_seqStats(
   double valDbl
){
   return ln_seqStats(valDbl) / 2.30258509299;
      /*convert ln(x) to log10(x)*/
      /*ln(10) = 2.30258509299*/
} /*log10_seqStats*/

/*-------------------------------------------------------\
| Fun03: pow10_seqStats
|   - rases 10 to the input exponant
|   - from: CRC standard mathmatical tables 26 edition
| Input:
|   - exponatDbl:
|     o rase base to power
| Output:
|   - Returns:
|     o 10^exponatDbl
\-------------------------------------------------------*/
double
pow10_seqStats(
   double exponatDbl
){
   /*x^a = e^[a * ln(x)]*/
   /*ln(x^r) = r * ln(x)*/
   /*e^x = 1 + x/1! + x^2/2! + x^3/3! + x^4/4! + ...*/

   double factoralDbl = 1;
   double retDbl = 1.0;
   double xPowerDbl = 0;
   signed int siEnd = 2;
   signed char negBl = 0;

   #define def_endDecimal_fun04 12 /*only do 10 rounds*/

   if(exponatDbl < 0)
   { /*If: negative power*/
      exponatDbl *= -1;
      negBl = 1;
   } /*If: negative power*/

   exponatDbl *= 2.30258509299; /*ln(10)*/

   /*find e^x*/
   retDbl = 1.0 + exponatDbl; /*1 + x/1!*/
   siEnd = 2;
   xPowerDbl = exponatDbl;
   while(siEnd < def_endDecimal_fun04)
   { /*Loop: e^x rest of series*/
      xPowerDbl *= exponatDbl;
      factoralDbl *= siEnd++;
      xPowerDbl /= factoralDbl;
      retDbl += xPowerDbl;
   } /*Loop: e^x rest of series*/

   if(negBl)
      return 1.0 / retDbl;
   else
      return retDbl;
} /*pow10_seqStats*/

/*-------------------------------------------------------\
| Fun04: qhistToMedQ_seqStats
|   - converts a q-score histogram to medain q-score
| Input:
|   - qhistArySI:
|     o signed int array with q-scores
|   - seqLenSI:
|     o length of sequence
|   - ontMethodBl:
|     o 1: use ONT method (convert Q to error, then back
|          to q-score with log)
|     o 0: use regular mean
| Output:
|   - Returns median q-score
\-------------------------------------------------------*/
float
qhistToMedQ_seqStats(
   signed int *qhistArySI,
   signed int seqLenSI,
   signed char ontMethodBl
){
   signed int siNt = 0;
   signed int qSI = 0;
   signed int midSI = 0;
   float medianQF = 0;

   midSI = seqLenSI >> 1;

   for(siNt = 0; siNt < def_maxQ_seqStats; ++siNt)
   { /*Loop: find median q-score*/
      qSI += qhistArySI[siNt];

      if(qSI >= midSI)
      { /*If: found the midpoint*/
         if(qSI > midSI || qSI & 1)
            medianQF = (float) siNt;
         else
         { /*Else: even & 2 differnt Q's at midpoint*/
            qSI = siNt++;
            while(qhistArySI[siNt++] == 0)
               ;
            --siNt; /*account for being one off*/
            medianQF = (qSI + siNt) / ((float) 2);
         } /*Else: even & 2 differnt Q's at midpoint*/

         break;
      } /*If: found the midpoint*/
   } /*Loop: find median q-score*/

   if(ontMethodBl)
   { /*If: using ont's q-score method*/
      medianQF /= (float) -10;

      #ifdef NO_MATH_H
         medianQF = pow10_seqStats(medianQF);
         medianQF = -10 * log10_seqStats(medianQF);
      #else
         medianQF = pow(10, medianQF);
         medianQF = -10 * log10(medianQF);
      #endif
   } /*If: using ont's q-score method*/

   return medianQF;
} /*qhistToMedQ_seqStats*/

/*-------------------------------------------------------\
| Fun05: qhistToMeanQ_seqStats
|   - converts a q-score histogram to mean q-score
| Input:
|   - qhistArySI:
|     o signed int array with q-scores
|   - seqLenSI:
|     o length of sequence
|   - ontMethodBl:
|     o 1: use ONT method (convert Q to error, then back
|          to q-score with log)
|     o 0: use regular mean
| Output:
|   - Returns mean q-score
\-------------------------------------------------------*/
float
qhistToMeanQ_seqStats(
   signed int *qhistArySI,
   signed int seqLenSI,
   signed char ontMethodBl
){
   signed int siNt = 0;
   double meanDbl = 0;
   float qF = 0;

   for(siNt = 0; siNt < def_maxQ_seqStats; ++siNt)
   { /*Loop: find median q-score*/
      if(! qhistArySI[siNt])
         ; /*no q-scores present*/

      else if(ontMethodBl)
      { /*Else If: using ont conversion*/
         qF = ((float) siNt / (float) -10);

         #ifdef NO_MATH_H
            qF = pow10_seqStats(qF);
         #else
            qF = pow(10, qF);
         #endif

         meanDbl += (qF * qhistArySI[siNt]) / seqLenSI;
      } /*Else If: using ont conversion*/

      else
         meanDbl += (qhistArySI[siNt] * siNt) / seqLenSI;
   } /*Loop: find median q-score*/

   if(ontMethodBl)
   { /*If: doing ont's conversion method*/
      #ifdef NO_MATH_H
         meanDbl = -10 * log10_seqStats(meanDbl);
      #else
         meanDbl = -10 * log10(meanDbl);
      #endif
   } /*If: doing ont's conversion method*/

   return meanDbl;
} /*qhistToMeanQ_seqStats*/

/*-------------------------------------------------------\
| Fun06: meanMedQ_seqStats
|   - gets the mean and median q-score for a sequence
| Input:
|   - qStr:
|     o c-string with q-score entry to get q-scores for
|   - seqLenSI:
|     o number of bases in the sequence (index 1)
|   - medianQFPtr:
|     o float pointer to get the median q-score
| Output:
|   - Modifies:
|     o medainQFPtr to have the median q-score
|   - Returns:
|     o mean q-score
\-------------------------------------------------------*/
float
meanMedQ_seqStats(
   signed char *qStr,
   signed int seqLenSI,
   float *medianQFPtr
){
   signed int histArySI[def_maxQ_seqStats];
   unsigned long totalQUL = 0;
   signed int siNt = 0;
   signed int qSI = 0;

   *medianQFPtr = 0;
   if(! qStr)
      return 0;
   else if(! qStr[0])
      return 0;

   for(siNt = 0; siNt < def_maxQ_seqStats; ++siNt)
      histArySI[siNt] = 0;

   for(siNt = 0; siNt < seqLenSI; ++siNt)
   { /*Loop: find the q-score*/
      qSI = qStr[siNt] - def_qAdj_seqStats;
      totalQUL += qSI;
      ++histArySI[qSI];
   } /*Loop: find the q-score*/

   *medianQFPtr =
      qhistToMedQ_seqStats(histArySI, seqLenSI, 0);
   return (double) totalQUL / (float) seqLenSI;
} /*meanMedQ_seqStats*/

/*-------------------------------------------------------\
| Fun07: ontMeanMedQ_seqStats
|   - gets the mean and median q-score for a sequence using
|     ONT's method
| Input:
|   - qStr:
|     o c-string with q-score entry to get q-scores for
|   - seqLenSI:
|     o number of bases in the sequence (index 1)
|   - medianQFPtr:
|     o float pointer to get the median q-score
|   - keepFirst60Bl:
|     o 1: keep the first 60 bases in the read
|     o 0: follow ONT and ingore first 60 bases
| Output:
|   - Modifies:
|     o medainQFPtr to have the median q-score
|   - Returns:
|     o mean q-score
\-------------------------------------------------------*/
float
ontMeanMedQ_seqStats(
   signed char *qStr,
   signed int seqLenSI,
   float *medianQFPtr,
   signed char keepFirst60Bl
){
   signed int histArySI[def_maxQ_seqStats];
   signed int siNt = 0;
   signed int qSI = 0;
   float qF = 0;

   *medianQFPtr = 0;
   if(! qStr)
      return 0;
   else if(! qStr[0])
      return 0;

   for(siNt = 0; siNt < def_maxQ_seqStats; ++siNt)
      histArySI[siNt] = 0;

   if(seqLenSI > 60 && ! keepFirst60Bl)
      siNt = 60;
   else
      siNt = 0;

   while(siNt < seqLenSI)
   { /*Loop: find the q-score*/
      qSI = qStr[siNt++] - def_qAdj_seqStats;
      ++histArySI[qSI];
   } /*Loop: find the q-score*/

   if(seqLenSI > 60 && ! keepFirst60Bl)
   { /*If: ignoring the first 60 bases*/
      *medianQFPtr =
         qhistToMedQ_seqStats(histArySI, seqLenSI-60, 1);
      qF = qhistToMeanQ_seqStats(histArySI,seqLenSI-60,1);
   } /*If: ignoring the first 60 bases*/

   else
   { /*Else: keeping the first 60 bases*/
      *medianQFPtr =
         qhistToMedQ_seqStats(histArySI, seqLenSI, 1);
            /*1 is for using ONT method*/
      qF = qhistToMeanQ_seqStats(histArySI, seqLenSI, 1);
   } /*Else: keeping the first 60 bases*/

   return qF;
} /*ontMeanMedQ_seqStats*/

/*-------------------------------------------------------\
| Fun08: pReadStats_seqStats
|   - prints stats for a single read
| Input:
|   - qStr:
|     o c-string with q-score entry to get q-scores for
|   - seqLenSI:
|     o number of bases in the sequence (index 1)
|   - idStr:
|     o c-string with id to print out
|   - pHeadBlPtr:
|     o 1: print header and set to 0
|     o 0: do not print header
|   - keepFirst60Bl:
|     o 1: for ONT mean keep the firt 60 bases
|     o 0: for ONT mean discard the first 60 bases
|   - statMedthodSC:
|     o 1: use ONT stat only
|     o 2: use mean/median (no error step) only
|     o 3: print both ONT and mean/median
|   - outFILE:
|     o FILE pionter to print stats to
| Output:
|   - Prints;
|     o read stats (and header if requested) to outFILE
|   - Modifies:
|     o pHeadBlPtr to be 0 (if is not 0)
\-------------------------------------------------------*/
void
pReadStats_seqStats(
   signed char *qStr,        /*quality score entry*/
   signed int seqLenSI,      /*bases in quality score*/
   signed char *idStr,       /*read id/name to print*/
   signed char *pHeadBlPtr,  /*1: print header*/
   signed char keepFirst60Bl,/*1: ONT keep first 60 nt*/
   signed char statMethodSC, /*1:ONT, 2:mean/med; 3:both*/
   void *outFILE             /*file to print to*/
){
   float meanF = 0;
   float medianF = 0;
   signed char outIdStr[126];
   signed int siNt = 0;


   if(*pHeadBlPtr)
   { /*If: printing the header*/
      fprintf((FILE *) outFILE, "id\tlength");
      if(statMethodSC & 1)
         fprintf(
            (FILE *) outFILE, 
            "\tont_mean_q\tont_median_q"
         );
      if(statMethodSC & 2)
         fprintf(
            (FILE *) outFILE, 
            "\tmean_q\tmedian_q%s",
            str_endLine
         );
      fprintf((FILE *) outFILE, "%s", str_endLine);
      *pHeadBlPtr = 0;
   } /*If: printing the header*/


   while(idStr[siNt] > 32)
   { /*Loop: copy read id*/
      outIdStr[siNt] = idStr[siNt];
      ++siNt;
   }  /*Loop: copy read id*/
   outIdStr[siNt] = 0;

   fprintf((FILE *) outFILE,"%s\t%i",outIdStr,seqLenSI);

   if(statMethodSC & 1)
   { /*If: printing ONT method*/
      meanF =
         ontMeanMedQ_seqStats(
            qStr,
            seqLenSI,
            &medianF,
            keepFirst60Bl
         );
      fprintf(
         (FILE *) outFILE,
         "\t%0.2f\t%0.2f",
         meanF,
         medianF
      );
   } /*If: printing ONT method*/

   if(statMethodSC & 2)
   { /*If: printing mean/median method*/
      meanF = meanMedQ_seqStats(qStr, seqLenSI, &medianF);
      fprintf(
         (FILE *) outFILE,
         "\t%0.2f\t%0.2f%s",
         meanF,
         medianF,
         str_endLine
      );
   } /*If: printing mean/median method*/

   else
      fprintf((FILE *) outFILE, "%s", str_endLine);
} /*pReadStats_seqStats*/
