/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' endLineBin SOF: Start Of File
'   - has line ending string for binaray's for different
'     OS's
'   - not really needed, but I want the flexability
'   - my code assumes there can be one or chararacters
'     for a break (no more)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef WINDOWS
   static signed char *str_endLineBin =
      (signed char *) "\r\n";
#endif

#ifdef OLD_MAC_BREAK
   static signed char *str_endLineBin =
      (signed char *) "\r";
#endif

#ifdef UNIX
   static signed char *str_endLineBin =
      (signed char *) "\n";
#endif

#ifdef LINUX
   static signed char *str_endLineBin =
      (signed char *) "\n";
#endif
