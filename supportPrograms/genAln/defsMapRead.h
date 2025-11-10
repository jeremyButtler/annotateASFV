/*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~\
' defsMapRead SOF: Start Of File
'   - default values for mapRead
'   o header:
'     - gaurds
'   o sec01:
'     - kmer scan variables
'   o sec02:
'     - chaining and scoring variables
\~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~*/

/*-------------------------------------------------------\
| Header:
|   - gaurds
\-------------------------------------------------------*/

#ifndef DEFAULTS_MAP_READ_H
#define DEFAULTS_MAP_READ_H

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec01:
^   - kmer scan variables
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

/*kmer scan defaults*/
#define def_minLen_defsMapRead 0
#define def_secLen_defsMapRead 100
#define def_thirdLen_defsMapRead 1000

#define def_minKmer_defsMapRead 5
#define def_secKmer_defsMapRead 10
#define def_thirdKmer_defsMapRead 13

#define def_minChain_defsMapRead 1
#define def_secChain_defsMapRead 1
#define def_thirdChain_defsMapRead 1


/* What these defaults mean:
` - my logic is keep is small for small references, then
`   expand it large for larger references
` - one: reference is < 100 bases
`   - 5mers; kmer space is 1024
`   - chains are one kmer long
` - two: reference is between 100 to 1000 bases
`   - 10mers, which should reduce matches, but might still
`     hit regions with high matches
`   - chains are at least 1 kmer (10 bases) long
` - third: reference is over 1000 bases
`   - 13mers (bit slower, but probably better job)
`   - chains are at least 1 kmer (13 bases) long
*/

/*>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>\
^ Sec02:
^   - chaining and scoring variables
\<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<*/

#define def_minScore_defsMapRead 0.4f
    /*min % score to keep alignment*/
#define def_minMatch_defsMapRead 0.0f
   /*at least 80% matches to keep read*/
#define def_minPercLen_defsMapRead 0.2f
    /*min % of query length mapped*/
#define def_chainMinLen_defsMapRead 0.20f
   /*min % query bases in merged chain*/
#define def_maxPercLen_defsMapRead 0.25f
    /*maximum gap length in percent query length to merge
    `  chains at
    */

#define def_alnEnd_defsMapRead 0
   /*1: is align softmasked ends
   `  gets stuff at end, but that is sometimes just
   `  nosie, so better to leave off
   */
#define def_subAln_defsMapRead 1 /*keep best sub-aln*/
#define def_chainGap_defsMapRead 0
  /*gap score for chaining step (0 is do not apply)*/
#define def_chainMatch_defsMapRead 5
  /*match score for chaining step*/

/*the Needleman/Waterman matrix's can be more complex,
`  however, in this case keep same as core.
`  Also, their scores are only used in gap filling,
`  not in the final score, see alnDefs for the defaults
*/

#endif
