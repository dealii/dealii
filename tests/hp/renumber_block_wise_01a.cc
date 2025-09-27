// ------------------------------------------------------------------------
//
// SPDX-License-Identifier: LGPL-2.1-or-later
// Copyright (C) 2009 - 2020 by the deal.II authors
//
// This file is part of the deal.II library.
//
// Part of the source code is dual licensed under Apache-2.0 WITH
// LLVM-exception OR LGPL-2.1-or-later. Detailed license information
// governing the source code and code contributions can be found in
// LICENSE.md and CONTRIBUTING.md at the top level directory of deal.II.
//
// ------------------------------------------------------------------------



// A redux of the _01 test that happened to fail on a branch to remove
// a bunch of iterator functions.



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

#include <deal.II/hp/fe_collection.h>

#include <deal.II/lac/vector.h>

#include "../tests.h"



template <int dim>
std::vector<types::global_dof_index>
get_dofs(const DoFHandler<dim> &dof)
{
  std::vector<types::global_dof_index> local;
  std::vector<types::global_dof_index> global;
  for (typename DoFHandler<dim>::active_cell_iterator cell = dof.begin_active();
       cell != dof.end();
       ++cell)
    {
      local.resize(cell->get_fe().dofs_per_cell);
      cell->get_dof_indices(local);

      global.insert(global.end(), local.begin(), local.end());
    }

  return global;
}



template <int dim>
void
check_renumbering(DoFHandler<dim> &dof)
{
  // do component-wise and save the
  // results
  DoFRenumbering::component_wise(dof);
  const std::vector<types::global_dof_index> vc = get_dofs(dof);

  deallog << "OK" << std::endl;
}


template <int dim>
void
check()
{
  Triangulation<dim> tr;
  GridGenerator::hyper_cube(tr, -1, 1);
  tr.refine_global(1);

  DoFHandler<dim> dof(tr);
  {
    bool coin = false;
    for (typename DoFHandler<dim>::active_cell_iterator cell =
           dof.begin_active();
         cell != dof.end();
         ++cell)
      {
        cell->set_active_fe_index(coin ? 0 : 1);
        coin = !coin;
      }
  }

  FESystem<dim> e1(FE_Q<dim>(2), 1, FE_DGQ<dim>(0), 1);
  FESystem<dim> e2(FE_Q<dim>(1), 1, FE_DGQ<dim>(0), 1);

  hp::FECollection<dim> fe_collection;
  fe_collection.push_back(e1);
  fe_collection.push_back(e2);

  dof.distribute_dofs(fe_collection);
  check_renumbering(dof);
  dof.clear();
}


int
main()
{
  initlog();
  deallog << std::setprecision(2);
  deallog << std::fixed;

  deallog.push("3d");
  check<3>();
  deallog.pop();
}
