// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// p4est defines a couple of function-like macros that are mixed in
// with other stuff in their header files so that in order to get to
// these macros, we would have to include the whole file. This of
// course defeats the purpose of wrapping things into module
// partitions.  Rather -- perhaps imprudently -- we repeat these
// macros here.

#ifndef dealii_p4est_macros_h
#define dealii_p4est_macros_h


// Taken from p4est.h, with a change from memset to std::memset:
#define P4EST_QUADRANT_INIT(q) \
  ((void)std::memset((q), -1, sizeof(p4est_quadrant_t)))
#define P8EST_QUADRANT_INIT(q) \
  ((void)std::memset((q), -1, sizeof(p4est_quadrant_t)))


#endif
