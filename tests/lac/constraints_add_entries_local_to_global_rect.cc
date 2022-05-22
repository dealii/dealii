// ---------------------------------------------------------------------
//
// Copyright (C) 2017 - 2022 by the deal.II authors
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



// this function tests the correctness of the implementation of
// AffineConstraints<double>::add_entries_local_to_global for row and column
// indices with different constraints on the rows and columns.

#include <deal.II/base/function.h>

#include <deal.II/lac/affine_constraints.h>
#include <deal.II/lac/dynamic_sparsity_pattern.h>

#include <iostream>

#include "../tests.h"

void
test()
{
  DynamicSparsityPattern    dsp1(4, 6), dsp2(4, 6);
  AffineConstraints<double> cm1, cm2;
  cm1.add_line(2);
  cm1.add_entry(2, 1, 0.5);
  cm1.add_entry(2, 3, 0.5);
  cm1.close();
  cm2.add_line(4);
  cm2.add_entry(4, 2, 0.3);
  cm2.add_entry(4, 5, 0.7);
  cm2.close();

  std::vector<types::global_dof_index> indices1(2), indices2(2);
  indices1[0] = 1;
  indices1[1] = 2;
  indices2[0] = 4;
  indices2[1] = 5;
  cm1.add_entries_local_to_global(indices1, indices2, dsp1);
  deallog << "Same constraints: " << std::endl;
  dsp1.print(deallog.get_file_stream());
  cm1.add_entries_local_to_global(indices1, cm2, indices2, dsp2);
  deallog << "Different constraints: " << std::endl;
  dsp2.print(deallog.get_file_stream());
}


int
main()
{
  initlog();

  test();
}
