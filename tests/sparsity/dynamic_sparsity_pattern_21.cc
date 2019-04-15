// ---------------------------------------------------------------------
//
// Copyright (C) 2004 - 2019 by the deal.II authors
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

// Check DynamicSparsityPattern::exists() with supplied rowindex by copying
// to a static SparsityPattern.

#include <deal.II/lac/dynamic_sparsity_pattern.h>
#include <deal.II/lac/sparsity_pattern.h>

#include "../tests.h"

int
main()
{
  initlog();

  const unsigned int n_dofs = 2;
  IndexSet           index_set(n_dofs);
  index_set.add_index(1);

  DynamicSparsityPattern dsp(n_dofs, n_dofs, index_set);
  SparsityPattern        sp;
  sp.copy_from(dsp);

  deallog << "OK" << std::endl;

  return 0;
}
