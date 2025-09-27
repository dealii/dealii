// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2022 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// check that conversion DoFHandler iterator to another DoFHandler iterator
// works


#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>

#include <deal.II/lac/full_matrix.h>

#include "../tests.h"

template <int dim, int spacedim>
void
test()
{
  Triangulation<dim, spacedim> triangulation;
  GridGenerator::hyper_cube(triangulation);
  triangulation.refine_global(2);
  FE_Q<dim, spacedim>       fe_1(2);
  FE_Q<dim, spacedim>       fe_2(3);
  DoFHandler<dim, spacedim> dof_handler_1(triangulation);
  DoFHandler<dim, spacedim> dof_handler_2(triangulation);
  dof_handler_1.distribute_dofs(fe_1);
  dof_handler_2.distribute_dofs(fe_2);

  std::vector<types::global_dof_index> local_dof_indices(
    fe_2.n_dofs_per_cell());

  for (const auto &it_dh_1 : dof_handler_1.active_cell_iterators())
    {
      const auto it_dh_2 = it_dh_1->as_dof_handler_iterator(dof_handler_2);
      Assert(it_dh_1->level() == it_dh_2->level(),
             ExcMessage("Iterator conversion failed: Level."));
      Assert(it_dh_1->index() == it_dh_2->index(),
             ExcMessage("Iterator conversion failed: Index."));
      Assert(it_dh_1->id() == it_dh_2->id(),
             ExcMessage("Iterator conversion failed: Id."));

      // Check that some basic features work (i.e. that we have the right
      // accessor type)
      it_dh_2->get_dof_indices(local_dof_indices);
    }
}



int
main()
{
  initlog();

  test<1, 1>();
  test<1, 2>();
  test<2, 2>();
  test<2, 3>();
  test<3, 3>();

  deallog << "OK" << std::endl;

  return 0;
}
