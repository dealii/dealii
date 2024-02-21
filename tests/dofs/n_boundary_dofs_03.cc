// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2017 - 2023 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// test that DoFHandler::n_boundary_dofs() yields the correct
// results in 1d, even if there are more than two boundary vertices
//
// same as the _02 test, but using the variant of the function that
// takes a std::set as argument


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

  // assign boundary ids
  triangulation.begin()->face(0)->set_boundary_id(12);
  triangulation.begin()->face(1)->set_boundary_id(13);
  std::next(triangulation.begin())->face(0)->set_boundary_id(14);
  std::next(triangulation.begin())->face(1)->set_boundary_id(15);


  FESystem<1, spacedim> fe(FE_Q<1, spacedim>(1), 1, FE_DGQ<1, spacedim>(1), 1);

  DoFHandler<1, spacedim> dof_handler(triangulation);

  dof_handler.distribute_dofs(fe);

  for (const types::boundary_id b : {12, 13, 14, 15})
    {
      const unsigned int N =
        dof_handler.n_boundary_dofs(std::set<types::boundary_id>{b});
      deallog << (int)b << ' ' << N << std::endl;
    }
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
