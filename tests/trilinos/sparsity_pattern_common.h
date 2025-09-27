// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check

#include <deal.II/lac/trilinos_sparsity_pattern.h>

#include "../tests.h"

#include "../testmatrix.h"

const unsigned int N = 15;


// reinitialize sparsity patterns for 5-point star
void
do_reinit(TrilinosWrappers::SparsityPattern &sp)
{
  sp.reinit((N - 1) * (N - 1), (N - 1) * (N - 1));
}



void
build_sparsity(TrilinosWrappers::SparsityPattern &sparsity_pattern)
{
  // generate usual 5-point sparsity pattern
  do_reinit(sparsity_pattern);
  FDMatrix(N, N).five_point_structure(sparsity_pattern);
  sparsity_pattern.compress();

  deallog << sparsity_pattern.n_rows() << " " << sparsity_pattern.n_cols()
          << " " << sparsity_pattern.bandwidth() << " "
          << sparsity_pattern.n_nonzero_elements() << std::endl;
}



void
row_length()
{
  TrilinosWrappers::SparsityPattern sparsity_pattern;
  build_sparsity(sparsity_pattern);

  for (unsigned int i = 0; i < sparsity_pattern.n_rows(); ++i)
    deallog << sparsity_pattern.row_length(i) << std::endl;

  deallog << "OK" << std::endl;
}



void
print_gnuplot()
{
  TrilinosWrappers::SparsityPattern sparsity_pattern;
  build_sparsity(sparsity_pattern);

  sparsity_pattern.print_gnuplot(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}



void
print()
{
  TrilinosWrappers::SparsityPattern sparsity_pattern;
  build_sparsity(sparsity_pattern);

  sparsity_pattern.print(deallog.get_file_stream());

  deallog << "OK" << std::endl;
}
