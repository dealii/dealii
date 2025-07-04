// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// since early 2009, the FEValues objects try to be more efficient by only
// recomputing things like gradients of shape functions if the cell on which
// we are is not a translation of the previous one. in this series of tests we
// make sure that this actually works the way it's supposed to be; in
// particular, if we create a mesh of two identical cells but one has a curved
// boundary, then they are the same if we use a Q1 mapping, but not a Q2
// mapping. so we test that the mass matrix we get from each of these cells is
// actually different in the latter case, but the same in the former
//
// Check cell to face lexicographic ordering


#include <deal.II/fe/fe_tools.h>

#include "../tests.h"



template <int dim>
void
test(unsigned int degree,
     unsigned int direction,
     bool         cell_hierarchical_numbering,
     bool         is_continuous)
{
  std::pair<std::vector<unsigned int>, std::vector<unsigned int>> result =
    FETools::cell_to_face_patch<dim>(degree,
                                     direction,
                                     cell_hierarchical_numbering,
                                     is_continuous);

  deallog << "Results for degree = " << degree << ", dim = " << dim
          << ", direction = " << direction
          << ", is_continuous = " << (is_continuous ? "true" : "false")
          << ", cell_hierarchical_numbering = "
          << (cell_hierarchical_numbering ? "true" : "false") << std::endl;

  deallog << "Cell 0: ";
  for (unsigned int val : result.first)
    deallog << val << " ";
  deallog << std::endl;

  deallog << "Cell 1: ";
  for (unsigned int val : result.second)
    deallog << val << " ";
  deallog << std::endl;
}


int
main()
{
  initlog();
  deallog << std::setprecision(3);

  test<2>(2, 0, false, false);
  test<2>(2, 0, false, true);
  test<2>(2, 0, true, true);

  test<2>(4, 0, false, true);

  test<2>(2, 1, false, false);
  test<2>(2, 1, false, true);


  test<3>(2, 0, false, false);
  test<3>(2, 0, false, true);
  test<3>(2, 0, true, true);

  test<3>(2, 1, false, false);
  test<3>(2, 1, false, true);

  test<3>(2, 2, false, false);
  test<3>(2, 2, false, true);
}
