// ---------------------------------------------------------------------
//
// Copyright (C) 1998 - 2013 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

#ifndef __deal2__magic_numbers_h
#define __deal2__magic_numbers_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

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


DEAL_II_NAMESPACE_CLOSE

#endif
