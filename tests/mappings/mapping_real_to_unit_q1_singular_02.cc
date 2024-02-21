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
// point which is close to a Cartesian mesh except for roundoff, which used to
// create divisions by zero in an earlier implementation

#include <deal.II/base/utilities.h>

#include <deal.II/fe/mapping_q1.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


void
test_real_to_unit_cell()
{
  const unsigned int dim = 2;

  Triangulation<dim> triangulation;

  std::vector<Point<dim>> points{
    Point<dim>(-0.29999999999999999, -0.29999999999999999),
    Point<dim>(-0.050000000000000003, -0.29999999999999999),
    Point<dim>(-0.29999999999999999, -0.050000000000000003),
    Point<dim>(-0.049999999999999989, -0.050000000000000003)};
  std::vector<CellData<dim>> cells(1);
  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].vertices[2] = 2;
  cells[0].vertices[3] = 3;
  cells[0].material_id = 0;
  triangulation.create_triangulation(points, cells, SubCellData());

  Point<dim>    point(-0.29999999999999999, -0.29999999999999999);
  MappingQ<dim> mapping(1);

  try
    {
      mapping.transform_real_to_unit_cell(triangulation.begin(), point);
    }
  catch (typename Mapping<dim>::ExcTransformationFailed &)
    {
      deallog << "Transformation for point " << point << " on cell with "
              << "center " << triangulation.begin()->center()
              << " is not invertible" << std::endl;
    }
  deallog << "OK" << std::endl;
}


int
main()
{
  initlog();

  test_real_to_unit_cell();

  return 0;
}
