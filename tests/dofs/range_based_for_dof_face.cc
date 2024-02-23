// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2018 - 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// Check face range-based for loops over DoFHandler objects

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>

#include <deal.II/fe/fe_q.h>

#include <deal.II/grid/grid_generator.h>

#include "../tests.h"


template <int dim, int spacedim>
void
check()
{
  Triangulation<dim, spacedim> tr;
  GridGenerator::hyper_cube(tr);
  tr.refine_global(2);
  DoFHandler<dim, spacedim> dh(tr);
  FE_Q<dim, spacedim>       fe(1);
  dh.distribute_dofs(fe);

  const unsigned int                   n_dofs_per_face = fe.n_dofs_per_face();
  std::vector<types::global_dof_index> face_dof_indices(n_dofs_per_face);

  for (const auto &cell : dh.active_cell_iterators())
    for (const auto &face : cell->face_iterators())
      {
        face->get_dof_indices(face_dof_indices);
        for (unsigned int i = 0; i < n_dofs_per_face; ++i)
          deallog << face_dof_indices[i] << ' ';
        deallog << std::endl;
      }
}


int
main()
{
  initlog();

  check<2, 2>();
  check<2, 3>();
  check<3, 3>();
}
