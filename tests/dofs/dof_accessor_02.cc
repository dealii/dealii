// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2013 - 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------

// Initialize a DOFAccessor for a vertex iterators and then read the
// dofs of the vertices.

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include "../tests.h"

int
main()
{
  initlog();
  Triangulation<2> triangulation;
  GridGenerator::hyper_cube(triangulation);
  DoFHandler<2>     dof_handler(triangulation);
  const FESystem<2> finite_element(FE_Q<2>(1), 2);
  dof_handler.distribute_dofs(finite_element);
  const auto vertex_end = triangulation.end_vertex();
  for (auto vertex = triangulation.begin_active_vertex(); vertex != vertex_end;
       ++vertex)
    {
      DoFAccessor<0, 2, 2, false> vertex_dofs(&triangulation,
                                              vertex->level(),
                                              vertex->index(),
                                              &dof_handler);
      deallog << vertex_dofs.dof_index(0) << ", " << vertex_dofs.dof_index(1)
              << std::endl;
    }
  return 0;
}
