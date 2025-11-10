/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' endLine SOF: Start Of File
'   - has line ending string for different OS's non-binary
'     fiels
'   - not really needed, but I want the flexability
'   - my code assumes there can be one or chararacters
'     for a break (no more)
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

#ifdef WIN_BREAK
   static signed char *str_endLine =
      (signed char *) "\r\n";
   #else
      #ifdef OLD_MAC_BREAK
         static signed char *str_endLine =
             (signed char *) "\r";
      #else
         #ifdef UNIX_BREAK
            static signed char *str_endLine =
               (signed char *) "\n";
         #else
            static signed char *str_endLine =
               (signed char *) "\n";
         #endif /*linux/default break*/
   #endif /*mac break*/
#endif /*windows break*/
