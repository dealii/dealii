// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2019 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// tests DoFRenumbering::random function for levels



#include <deal.II/base/function_lib.h>

#include <deal.II/dofs/dof_accessor.h>
#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>

#include <deal.II/fe/fe_dgp.h>
#include <deal.II/fe/fe_dgq.h>
#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>

#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/tria.h>
#include <deal.II/grid/tria_iterator.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
void
print_dofs(const DoFHandler<dim> &dof, unsigned int level)
{
  std::vector<types::global_dof_index> v(dof.get_fe().dofs_per_cell);
  for (typename DoFHandler<dim>::cell_iterator cell = dof.begin(level);
       cell != dof.end(level);
       ++cell)
    {
      deallog << "Cell " << cell << " -- ";
      cell->get_mg_dof_indices(v);
      for (unsigned int i = 0; i < v.size(); ++i)
        deallog << v[i] << ' ';
      deallog << std::endl;
    }
}


template <int dim>
void
test()
{
  Triangulation<dim> tr(Triangulation<dim>::limit_level_difference_at_vertices);
  GridGenerator::hyper_cube(tr);
  tr.refine_global((dim == 1 ? 2 : 1));

  DoFHandler<dim> dof_handler(tr);

  FE_Q<dim> fe(1);
  dof_handler.distribute_dofs(fe);
  dof_handler.distribute_mg_dofs();


  // Print dofs before reordering
  deallog << "Before reorder: " << std::endl;
  for (unsigned int level = 0; level < tr.n_levels(); ++level)
    {
      deallog << "Level " << level << std::endl;
      print_dofs(dof_handler, level);
    }

  for (unsigned int level = 0; level < tr.n_levels(); ++level)
    DoFRenumbering::random(dof_handler, level);

  // Print dofs after reordering
  deallog << "After reorder: " << std::endl;
  for (unsigned int level = 0; level < tr.n_levels(); ++level)
    {
      deallog << "Level " << level << std::endl;
      print_dofs(dof_handler, level);
    }
}


int
main()
{
  initlog();

  deallog << "1D" << std::endl;
  test<1>();
  deallog << std::endl << "2D" << std::endl;
  test<2>();
  deallog << std::endl << "3D" << std::endl;
  test<3>();
}
