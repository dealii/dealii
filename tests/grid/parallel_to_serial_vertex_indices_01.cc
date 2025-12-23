// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2015 - 2025 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------


// GridTools::parallel_to_serial_vertex_indices() test


#include <deal.II/distributed/fully_distributed_tria.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"


template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> tria;
  GridGenerator::subdivided_hyper_cube(tria, 5 - dim);

  parallel::fullydistributed::Triangulation<dim, spacedim> distributed_tria;
  distributed_tria.copy_triangulation(tria);

  const auto vertex_indices =
    GridTools::parallel_to_serial_vertex_indices(tria, distributed_tria);

  unsigned int vertex_id = 0;
  for (const auto &v = distributed_tria.get_vertices())
    {
    }



  int main()
  {
    initlog();

    test();

    return 0;
  }
