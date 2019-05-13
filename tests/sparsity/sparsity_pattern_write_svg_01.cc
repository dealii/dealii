// ---------------------------------------------------------------------
//
// Copyright (C) 2008 - 2018 by the deal.II authors
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
