// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

/*
 * Workaround for a bug in the Xcode generator.
 *
 * This file contains a dummy global symbol to trigger the link phase in
 * the generated Xcode project.
 */

const int global_symbol_42{42};
void
use_global_symbol_42()
{
  (void)global_symbol_42;
}

/**
 * If we are running the contrib/utilities/run_clang_tidy.sh script, we
 * generate a header file allheaders.h that includes all deal.II
 * headers. Include the file here so that all headers are being checked.
 */
#ifdef CLANG_TIDY
#  include <deal.II/allheaders.h>
#endif
