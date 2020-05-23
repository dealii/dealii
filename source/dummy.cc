// ---------------------------------------------------------------------
//
// Copyright (C) 2013 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE.md at
// the top level directory of deal.II.
//
// ---------------------------------------------------------------------

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
