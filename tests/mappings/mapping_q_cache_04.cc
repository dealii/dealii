// ---------------------------------------------------------------------
//
// Copyright (C) 2019 by the deal.II authors
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

// Test MappingQCache initialization with lambda

#include <deal.II/fe/fe_nothing.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_values.h>
#include <deal.II/fe/mapping_q_cache.h>
#include <deal.II/fe/mapping_q_generic.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

template <int dim>
void
do_test(const unsigned int degree)
{
  Triangulation<dim> tria;
  if (dim > 1)
    GridGenerator::hyper_ball(tria);
  else
    GridGenerator::hyper_cube(tria, -1, 1);

  MappingQGeneric<dim> mapping(degree);
  MappingQCache<dim>   mapping_cache(degree);
  Point<dim>           shift;
  for (unsigned int d = 0; d < dim; ++d)
    shift[d] = -0.5 + 0.1 * d;

  FE_Q<dim> fe_q(degree);

  const auto position_lambda =
    [&](const typename Triangulation<dim>::cell_iterator &cell)
    -> std::vector<Point<dim>> {
    FE_Nothing<dim> fe;
    FEValues<dim>   fe_values(mapping,
                            fe,
                            Quadrature<dim>(fe_q.get_unit_support_points()),
                            update_quadrature_points);

    std::vector<Point<dim>> support_points_moved(fe_q.dofs_per_cell);
    fe_values.reinit(cell);
    for (unsigned int i = 0; i < fe_q.dofs_per_cell; ++i)
      {
        const Point<dim> point  = fe_values.quadrature_point(i);
        Point<dim>       out    = point + shift;
        support_points_moved[i] = out;
      }
    return support_points_moved;
  };

  mapping_cache.initialize(tria, position_lambda);

  Point<dim> p1;
  for (unsigned int d = 0; d < dim; ++d)
    p1[d] = 0.2 + d * 0.15;
  Point<dim> p2;
  for (unsigned int d = 0; d < dim; ++d)
    p2[d] = 0.5;

  deallog << "Testing degree " << degree << " in " << dim << "D" << std::endl;
  for (const auto &cell : tria.cell_iterators())
    {
      deallog << "cell " << cell->id() << ": "
              << mapping_cache.transform_unit_to_real_cell(cell, p1)
              << " vs reference "
              << mapping.transform_unit_to_real_cell(cell, p1) + shift
              << std::endl;
      deallog << "cell " << cell->id() << ": "
              << mapping_cache.transform_unit_to_real_cell(cell, p2)
              << " vs reference "
              << mapping.transform_unit_to_real_cell(cell, p2) + shift
              << std::endl;
      const Point<dim> p_back = mapping_cache.transform_real_to_unit_cell(
        cell, mapping_cache.transform_unit_to_real_cell(cell, p2));
      deallog << "cell " << cell->id() << " pull-back: " << p_back << std::endl;
      AssertThrow((p2 - p_back).norm() < 1e-13, ExcInternalError());
    }
  deallog << std::endl;

  tria.refine_global(1);
  mapping_cache.initialize(tria, position_lambda);

  deallog << "Testing degree " << degree << " in " << dim << "D" << std::endl;
  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << "cell " << cell->id() << ": "
              << mapping_cache.transform_unit_to_real_cell(cell, p1)
              << " vs reference "
              << mapping.transform_unit_to_real_cell(cell, p1) + shift
              << std::endl;
      deallog << "cell " << cell->id() << ": "
              << mapping_cache.transform_unit_to_real_cell(cell, p2)
              << " vs reference "
              << mapping.transform_unit_to_real_cell(cell, p2) + shift
              << std::endl;
      const Point<dim> p_back = mapping_cache.transform_real_to_unit_cell(
        cell, mapping_cache.transform_unit_to_real_cell(cell, p2));
      deallog << "cell " << cell->id() << " pull-back: " << p_back << std::endl;
      AssertThrow((p2 - p_back).norm() < 1e-13, ExcInternalError());
    }
  deallog << std::endl;
}


int
main()
{
  initlog();
  do_test<1>(1);
  do_test<1>(3);
  do_test<2>(1);
  do_test<2>(3);
  do_test<2>(4);
  do_test<3>(1);
  do_test<3>(2);
  do_test<3>(3);
}
