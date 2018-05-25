// ---------------------------------------------------------------------
//
// Copyright (C) 2017 by the deal.II authors
//
// This file is part of the deal.II library.
//
// The deal.II library is free software; you can use it, redistribute
// it, and/or modify it under the terms of the GNU Lesser General
// Public License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
// The full text of the license can be found in the file LICENSE at
// the top level of the deal.II distribution.
//
// ---------------------------------------------------------------------

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
