/*----------------------------   magic_numbers.h     ---------------------------*/
/*      $Id$                 */
#ifndef __magic_numbers_H
#define __magic_numbers_H
/*----------------------------   magic_numbers.h     ---------------------------*/


// This is a list of magic numbers used throughout the library.
// They are collected in one file to avoid double usage.
// Naming convention: all names have to start with the sequence
// "mn_" denoting a magic number, then the library part follows
// (e.g. "tria_" or "dof_") and finally the purpose.

const unsigned int mn_tria_refine_flags_begin    = 0xabcc;
const unsigned int mn_tria_refine_flags_end      = 0xabcd;
const unsigned int mn_tria_line_user_flags_begin = 0xabce;
const unsigned int mn_tria_line_user_flags_end   = 0xabcf;
const unsigned int mn_tria_quad_user_flags_begin = 0xabce;
const unsigned int mn_tria_quad_user_flags_end   = 0xabcf;



/*----------------------------   magic_numbers.h     ---------------------------*/
/* end of #ifndef __magic_numbers_H */
#endif
/*----------------------------   magic_numbers.h     ---------------------------*/
