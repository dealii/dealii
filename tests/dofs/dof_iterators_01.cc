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



// check that conversion between Triangulation iterator and DoFHandler iterator
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
  FE_Q<dim, spacedim>       fe(2);
  DoFHandler<dim, spacedim> dof_handler(triangulation);
  dof_handler.distribute_dofs(fe);

  std::vector<types::global_dof_index> local_dof_indices(fe.n_dofs_per_cell());

  for (const auto &it_tria : triangulation.active_cell_iterators())
    {
      const auto it_dh = it_tria->as_dof_handler_iterator(dof_handler);
      Assert(it_tria->level() == it_dh->level(),
             ExcMessage("Iterator conversion failed: Level."));
      Assert(it_tria->index() == it_dh->index(),
             ExcMessage("Iterator conversion failed: Index."));
      Assert(it_tria->id() == it_dh->id(),
             ExcMessage("Iterator conversion failed: Id."));

      // Check that some basic features work (i.e. that we have the right
      // accessor type)
      it_dh->get_dof_indices(local_dof_indices);
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
