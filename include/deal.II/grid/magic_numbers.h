// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 1998 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#ifndef dealii_magic_numbers_h
#define dealii_magic_numbers_h

#include <deal.II/base/config.h>

DEAL_II_NAMESPACE_OPEN

// This is a list of magic numbers used throughout the library.
// They are collected in one file to avoid double usage.
// Naming convention: all names have to start with the sequence
// "mn_" denoting a magic number, then the library part follows
// (e.g. "tria_" or "dof_") and finally the purpose.

const unsigned int mn_tria_refine_flags_begin     = 0xa000;
const unsigned int mn_tria_refine_flags_end       = 0xa001;
const unsigned int mn_tria_coarsen_flags_begin    = 0xa010;
const unsigned int mn_tria_coarsen_flags_end      = 0xa011;
const unsigned int mn_tria_line_user_flags_begin  = 0xa100;
const unsigned int mn_tria_line_user_flags_end    = 0xa101;
const unsigned int mn_tria_quad_user_flags_begin  = 0xa110;
const unsigned int mn_tria_quad_user_flags_end    = 0xa111;
const unsigned int mn_tria_hex_user_flags_begin   = 0xa112;
const unsigned int mn_tria_hex_user_flags_end     = 0xa113;
const unsigned int mn_persistent_tria_flags_begin = 0xa200;
const unsigned int mn_persistent_tria_flags_end   = 0xa201;


DEAL_II_NAMESPACE_CLOSE

#endif
