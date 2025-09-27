// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

#include <deal.II/multigrid/sparse_matrix_collection.h>

#include <algorithm>

#include "../tests.h"



int
main()
{
  initlog();

  mg::SparseMatrixCollection<float> smc;
  smc.resize(0, 5);

  deallog << "matrix      " << smc.matrix.min_level() << '-'
          << smc.matrix.max_level() << std::endl;
  deallog << "matrix_in   " << smc.matrix_in.min_level() << '-'
          << smc.matrix_in.max_level() << std::endl;
  deallog << "matrix_out  " << smc.matrix_out.min_level() << '-'
          << smc.matrix_out.max_level() << std::endl;
  deallog << "matrix_up   " << smc.matrix_up.min_level() << '-'
          << smc.matrix_up.max_level() << std::endl;
  deallog << "matrix_down " << smc.matrix_down.min_level() << '-'
          << smc.matrix_down.max_level() << std::endl;
}
