/*----------------------------   magic_numbers.h     ---------------------------*/
/*      $Id$                 */
/*      Copyright W. Bangerth, University of Heidelberg, 1998 */
#ifndef __magic_numbers_H
#define __magic_numbers_H
/*----------------------------   magic_numbers.h     ---------------------------*/


// This is a list of magic numbers used throughout the library.
// They are collected in one file to avoid double usage.
// Naming convention: all names have to start with the sequence
// "mn_" denoting a magic number, then the library part follows
// (e.g. "tria_" or "dof_") and finally the purpose.

const unsigned int mn_tria_refine_flags_begin    = 0xa000;
const unsigned int mn_tria_refine_flags_end      = 0xa001;
const unsigned int mn_tria_coarsen_flags_begin   = 0xa010;
const unsigned int mn_tria_coarsen_flags_end     = 0xa011;
const unsigned int mn_tria_line_user_flags_begin = 0xa100;
const unsigned int mn_tria_line_user_flags_end   = 0xa101;
const unsigned int mn_tria_quad_user_flags_begin = 0xa110;
const unsigned int mn_tria_quad_user_flags_end   = 0xa111;
const unsigned int mn_tria_hex_user_flags_begin  = 0xa112;
const unsigned int mn_tria_hex_user_flags_end    = 0xa113;



/*----------------------------   magic_numbers.h     ---------------------------*/
/* end of #ifndef __magic_numbers_H */
#endif
/*----------------------------   magic_numbers.h     ---------------------------*/
