// ---------------------------------------------------------------------
//
// Copyright (C) 2014 - 2020 by the deal.II authors
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


// Test by Nicola Giuliani: compute the measure of a 2d distorted cell in codim
// one

#include <deal.II/base/quadrature_lib.h>

#include <deal.II/dofs/dof_accessor.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_out.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include "../tests.h"



void
test()
{
  Triangulation<2, 3>      tria;
  std::vector<Point<3>>    vertices;
  std::vector<CellData<2>> cells;
  SubCellData              subcelldata;

  const double tol = 1e-12;
  vertices.push_back(Point<3>{0, 0, 1});
  vertices.push_back(Point<3>{1, 0, -10});
  vertices.push_back(Point<3>{0, 1, -1});
  vertices.push_back(Point<3>{1, 1, 1});

  cells.resize(1);

  cells[0].vertices[0] = 0;
  cells[0].vertices[1] = 1;
  cells[0].vertices[2] = 2;
  cells[0].vertices[3] = 3;

  tria.create_triangulation(vertices, cells, subcelldata);

  FE_Q<2, 3>       fe(1);
  DoFHandler<2, 3> dof_handler(tria);
  dof_handler.distribute_dofs(fe);

  QGauss<2>       quadrature_formula(4);
  MappingQ1<2, 3> mapping;
  FEValues<2, 3>  fe_values(mapping, fe, quadrature_formula, update_JxW_values);

  double sum_2 = 0;
  auto   cell  = dof_handler.begin_active();

  const double sum_1 = cell->measure();
  fe_values.reinit(cell);
  for (unsigned int q = 0; q < quadrature_formula.size(); ++q)
    {
      sum_2 += fe_values.JxW(q);
    }
  if (std::abs(sum_1 - sum_2) < tol)
    deallog << "Test passed" << std::endl;
  else
    deallog << sum_1 << " " << sum_2 << std::endl;
  deallog << std::endl;
}



int
main()
{
  initlog();
  deallog << std::setprecision(8);


  test();
}
