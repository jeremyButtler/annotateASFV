){
   signed int adjustSI = 0;
   double retDbl = 0;
   double squaredDbl = 0;
   signed int roundSI = 0;

   #define EULER_CONST_fun03 2.718281828459045235f
   #define TAYLOR_ITERATIONS_fun03 20
   #define log_10_of_2_fun03 0.30102999566f

   if(valDbl <= 0)
      return 0.0/0.0;

   /*get whole number part out*/
   while(valDbl >= 10)
   { /*Loop: reduce variable to no powers of 10*/
      valDbl /= 10;
      ++adjustSI;
   } /*Loop: reduce variable to no powers of 10*/

   while(valDbl < 1)
   { /*Loop: move variable to less then 10*/
      valDbl *= 10;
      --adjustSI;
   } /*Loop: move variable to less then 10*/

   /*find natural log*/
   valDbl = (valDbl - 1) * (valDbl + 1);
   retDbl = adjustSI;
   squaredDbl = valDbl * valDbl;
   roundSI = 1;

   while(roundSI < TAYLOR_ITERATIONS_fun03)
   { /*Loop: find log2*/
      retDbl += ( valDbl / ((roundSI + 1) << 1) );
      valDbl *= squaredDbl;
      ++roundSI;
   }  /*Loop: find log2*/

   return retDbl;
} /*log10_seqStats*/
