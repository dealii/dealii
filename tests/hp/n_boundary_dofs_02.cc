// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2010 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that hp::DoFHandler::n_boundary_dofs() yields the correct
// results in 1d, even if there are more than two boundary vertices


#include <deal.II/base/tensor.h>
#include <deal.II/base/utilities.h>

#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/intergrid_map.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_accessor.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/hp/fe_collection.h>

#include "../tests.h"



template <int spacedim>
void
test()
{
  // create a triangulation that spans the disjoint interval [0,1] \union [2,3]
  Triangulation<1, spacedim> triangulation_1;
  GridGenerator::hyper_cube(triangulation_1, 0, 1);

  Triangulation<1, spacedim> triangulation_2;
  GridGenerator::hyper_cube(triangulation_2, 2, 3);

  Triangulation<1, spacedim> triangulation;
  GridGenerator::merge_triangulations(triangulation_1,
                                      triangulation_2,
                                      triangulation);


  hp::FECollection<1, spacedim> fe;
  fe.push_back(
    FESystem<1, spacedim>(FE_Q<1, spacedim>(1), 1, FE_DGQ<1, spacedim>(1), 1));
  fe.push_back(FESystem<1, spacedim>(FE_Q<1, spacedim>(2), 2));

  DoFHandler<1, spacedim> dof_handler(triangulation);

  unsigned int index = 0;
  for (typename DoFHandler<1, spacedim>::active_cell_iterator cell =
         dof_handler.begin_active();
       cell != dof_handler.end();
       ++cell, ++index)
    cell->set_active_fe_index(index);

  dof_handler.distribute_dofs(fe);

  const unsigned int N = dof_handler.n_boundary_dofs();
  deallog << N << std::endl;
}


int
main()
{
  initlog();

  deallog.push("2d");
  test<2>();
  deallog.pop();

  deallog.push("3d");
  test<3>();
  deallog.pop();
}
