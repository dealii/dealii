// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Test MappingQCache by comparison with MappingQ in parallel

#include <deal.II/base/mpi.h>

#include <deal.II/distributed/tria.h>

#include <deal.II/fe/mapping_q.h>
#include <deal.II/fe/mapping_q_cache.h>

#include <deal.II/grid/grid_generator.h>

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

  tria.refine_global(1);

  MappingQ<dim>      mapping(degree);
  MappingQCache<dim> mapping_cache(degree);
  mapping_cache.initialize(mapping, tria);

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
}


int
main(int argc, char **argv)
{
  Utilities::MPI::MPI_InitFinalize mpi_initialization(argc, argv, 1);
  MPILogInitAll                    log;

  do_test<1>(1);
  do_test<1>(3);
  do_test<2>(1);
  do_test<2>(3);
  do_test<2>(4);
  do_test<3>(1);
  do_test<3>(2);
  do_test<3>(3);
}
