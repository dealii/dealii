// ---------------------------------------------------------------------
//
// Copyright (C) 2006 - 2017 by the deal.II authors
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

// check that MappingQ1::transform_real_to_unit_cell can handle the case of a
// singular discriminant. Previously we used to have a division by zero

#include "../tests.h"

#include <deal.II/base/utilities.h>
#include <deal.II/fe/mapping_q1.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

template <int dim>
void
test_real_to_unit_cell()
{
  deallog << "dim=" << dim << std::endl;

  Triangulation<dim> triangulation;
  GridGenerator::hyper_ball(triangulation);

  Point<dim>           point;
  MappingQGeneric<dim> mapping(1);

  point[1] = -1. / (1 + std::sqrt(2.0)) / std::sqrt(2);

  // check on cell 2
  typename Triangulation<dim>::cell_iterator cell = triangulation.begin();
  ++cell;
  try
    {
      mapping.transform_real_to_unit_cell(cell, point);
    }
  catch(typename Mapping<dim>::ExcTransformationFailed)
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
  std::ofstream logfile("output");
  deallog.attach(logfile);

  test_real_to_unit_cell<2>();

  return 0;
}
