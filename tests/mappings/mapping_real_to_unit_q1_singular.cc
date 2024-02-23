// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2021 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that MappingQ1::transform_real_to_unit_cell can handle the case of a
// singular discriminant. Previously we used to have a division by zero

#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim>
void
test_real_to_unit_cell()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball(triangulation);

  Point<dim>    point;
  MappingQ<dim> mapping(1);

  point[1] = -1. / (1 + std::sqrt(2.0)) / std::sqrt(2);

  // check on cell 2
  typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
  ++cell;
  try
    {
      mapping.transform_real_to_unit_cell(cell, point);
    }
  catch (const typename Mapping<dim>::ExcTransformationFailed &)
    {
      deallog << "Transformation for point " << point << " on cell with "
              << "center " << cell->center() << " is not invertible"
              << std::endl;
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test_real_to_unit_cell<2>();
}
