// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check SparsityPattern::write_svg

#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

int
main()
{
  initlog();

  SparsityPattern sparsity(2, 3, 3);
  for (unsigned int i = 0; i < 2; ++i)
    for (unsigned int j = 0; j < 3; ++j)
      if (i < j)
        sparsity.add(i, j);
  sparsity.compress();

  sparsity.print_svg(deallog.get_file_stream());
}
