//----------------------------  magic_numbers.h  ---------------------------
//    $Id$
//    Version: $Name$
//
//    Copyright (C) 1998, 1999, 2000, 2001, 2002 by the deal.II authors
//
//    This file is subject to QPL and may not be  distributed
//    without copyright and license information. Please refer
//    to the file deal.II/doc/license.html for the  text  and
//    further information on this license.
//
//----------------------------  magic_numbers.h  ---------------------------
#ifndef __deal2__magic_numbers_h
#define __deal2__magic_numbers_h

#include <base/config.h>


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
const unsigned int mn_persistent_tria_flags_begin= 0xa200;
const unsigned int mn_persistent_tria_flags_end  = 0xa201;


#endif
