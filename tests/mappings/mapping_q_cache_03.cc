// ---------------------------------------------------------------------
//
// Copyright (C) 2019 - 2022 by the deal.II authors
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

// Test MappingQCache by comparison with MappingQ for the case when we
// change the mesh

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_cache.h>

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

  MappingQ<dim>      mapping(degree);
  MappingQCache<dim> mapping_cache(degree);
  mapping_cache.initialize(tria, mapping);

  Point<dim> p1;
  for (unsigned int d = 0; d < dim; ++d)
    p1[d] = 0.2 + d * 0.15;
  Point<dim> p2;
  for (unsigned int d = 0; d < dim; ++d)
    p2[d] = 0.5;

  deallog << "Testing degree " << degree << " in " << dim << 'D' << std::endl;
  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << "cell " << cell->id() << ": "
              << mapping_cache.transform_unit_to_real_cell(cell, p1)
              << " vs reference "
              << mapping.transform_unit_to_real_cell(cell, p1) << std::endl;
      deallog << "cell " << cell->id() << ": "
              << mapping_cache.transform_unit_to_real_cell(cell, p2)
              << " vs reference "
              << mapping.transform_unit_to_real_cell(cell, p2) << std::endl;
      AssertThrow((p2 -
                   mapping_cache.transform_real_to_unit_cell(
                     cell, mapping_cache.transform_unit_to_real_cell(cell, p2)))
                      .norm() < 1e-13,
                  ExcInternalError());
    }
  deallog << std::endl;

  tria.refine_global(1);
  mapping_cache.initialize(tria, mapping);

  deallog << "Testing degree " << degree << " in " << dim << 'D' << std::endl;
  for (const auto &cell : tria.active_cell_iterators())
    {
      deallog << "cell " << cell->id() << ": "
              << mapping_cache.transform_unit_to_real_cell(cell, p1)
              << " vs reference "
              << mapping.transform_unit_to_real_cell(cell, p1) << std::endl;
      deallog << "cell " << cell->id() << ": "
              << mapping_cache.transform_unit_to_real_cell(cell, p2)
              << " vs reference "
              << mapping.transform_unit_to_real_cell(cell, p2) << std::endl;
      AssertThrow((p2 -
                   mapping_cache.transform_real_to_unit_cell(
                     cell, mapping_cache.transform_unit_to_real_cell(cell, p2)))
                      .norm() < 1e-13,
                  ExcInternalError());
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
