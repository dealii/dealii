// -----------------------------------------------------------------------------
//
// SPDX-License-Identifier: Apache-2.0 WITH LLVM-exception OR LGPL-2.1-or-later
// Copyright (C) 2017 - 2018 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Detailed license information governing the source code and contributions
// can be found in LICENSE.md and CONTRIBUTING.md at the top level directory.
//
// -----------------------------------------------------------------------------

// Check DoFTools::write_gnuplot_dof_support_point_info for a particular set of
// points

#include <deal.II/dofs/dof_tools.h>

#include "../tests.h"

void
test(const double epsilon)
{
  constexpr unsigned int dim = 2;
  Point<dim>             p1(1.0, 1.0);
  Point<dim>             p2(1.0 + epsilon, 0.0);

  std::map<types::global_dof_index, Point<dim>> support_points;
  support_points[0] = p1;
  support_points[1] = p2;

  DoFTools::write_gnuplot_dof_support_point_info(deallog.get_file_stream(),
                                                 support_points);
  deallog << std::endl;
}



int
main()
{
  initlog();
  test(1e-4);
  test(1e-5);
  test(1e-6);
}
